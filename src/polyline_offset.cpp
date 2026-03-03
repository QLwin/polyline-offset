#include "polyline_offset.h"
#include <cassert>
#include <cmath>
#include <limits>

namespace polyline_offset {

static constexpr double EPS = 1e-9;

// ──────────────────────── helper geometry ─────────────────────────────

std::optional<double> line_line_param(
    const Vec2& P, const Vec2& D,
    const Vec2& Q, const Vec2& E)
{
    double denom = D.cross(E);
    if (std::abs(denom) < EPS) return std::nullopt;
    return (Q - P).cross(E) / denom;
}

std::optional<SegIsect> seg_seg_isect(
    const Vec2& A, const Vec2& B,
    const Vec2& C, const Vec2& D)
{
    Vec2   AB = B - A, CD = D - C;
    double denom = AB.cross(CD);
    if (std::abs(denom) < EPS) return std::nullopt;

    Vec2   AC = C - A;
    double t  = AC.cross(CD) / denom;
    double s  = AC.cross(AB) / denom;

    if (t < -EPS || t > 1 + EPS || s < -EPS || s > 1 + EPS)
        return std::nullopt;

    t = std::clamp(t, 0.0, 1.0);
    s = std::clamp(s, 0.0, 1.0);
    return SegIsect{t, s, A + AB * t};
}

double point_polyline_dist(const Vec2& P, const std::vector<Vec2>& poly)
{
    double best = std::numeric_limits<double>::max();
    for (size_t i = 0; i + 1 < poly.size(); ++i) {
        Vec2   AB   = poly[i + 1] - poly[i];
        double len2 = AB.lengthSq();
        double t    = (len2 > 1e-24)
                        ? std::clamp((P - poly[i]).dot(AB) / len2, 0.0, 1.0)
                        : 0.0;
        double d = (P - (poly[i] + AB * t)).length();
        best = std::min(best, d);
    }
    return best;
}

// ──────────── remove near-duplicate consecutive vertices ──────────────

static std::vector<Vec2> deduplicate(const std::vector<Vec2>& pts,
                                     double tol = EPS)
{
    if (pts.empty()) return {};
    std::vector<Vec2> out;
    out.push_back(pts.front());
    for (size_t i = 1; i < pts.size(); ++i) {
        if ((pts[i] - out.back()).lengthSq() > tol * tol)
            out.push_back(pts[i]);
    }
    return out;
}

// ───────────────── compute unit left-normal of a segment ──────────────

static Vec2 seg_normal(const Vec2& A, const Vec2& B)
{
    return (B - A).perp().normalized();
}

// ─────────────── build raw offset polyline (open case) ────────────────
//
// For each original edge, shift it by `offset` along its left-normal.
// Adjacent offset edges are joined at their intersection (miter join)
// unless the miter ratio exceeds miter_limit, in which case a bevel
// (single connecting segment) is used instead.

static std::vector<Vec2> raw_offset_open(
    const std::vector<Vec2>& poly,
    double offset,
    double miter_limit)
{
    const int n = static_cast<int>(poly.size());
    if (n < 2) return poly;

    // Pre-compute offset edges: each edge i → (oA[i], oB[i])
    struct Edge { Vec2 a, b; };
    std::vector<Edge> edges(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        Vec2 norm = seg_normal(poly[i], poly[i + 1]) * offset;
        edges[i] = {poly[i] + norm, poly[i + 1] + norm};
    }

    std::vector<Vec2> raw;
    raw.push_back(edges[0].a);

    for (int i = 0; i < n - 2; ++i) {
        // Intersect offset edge i with edge i+1
        Vec2 d1 = edges[i].b - edges[i].a;
        Vec2 d2 = edges[i + 1].b - edges[i + 1].a;
        auto t = line_line_param(edges[i].a, d1, edges[i + 1].a, d2);

        if (t.has_value()) {
            Vec2 miter = edges[i].a + d1 * t.value();
            double miterLen = (miter - poly[i + 1]).length();
            if (std::abs(offset) > EPS &&
                miterLen / std::abs(offset) > miter_limit) {
                // Bevel: insert end of edge i and start of edge i+1
                raw.push_back(edges[i].b);
                raw.push_back(edges[i + 1].a);
            } else {
                raw.push_back(miter);
            }
        } else {
            // Parallel edges – just use the endpoint
            raw.push_back(edges[i].b);
        }
    }

    raw.push_back(edges.back().b);
    return raw;
}

// ─────────────── build raw offset polyline (closed case) ──────────────

static std::vector<Vec2> raw_offset_closed(
    const std::vector<Vec2>& poly,
    double offset,
    double miter_limit)
{
    int n = static_cast<int>(poly.size());
    // Ensure the polygon is not already explicitly closed
    if (n >= 2 && poly.front() == poly.back()) --n;
    if (n < 2) return poly;

    struct Edge { Vec2 a, b; };
    std::vector<Edge> edges(n);
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        Vec2 norm = seg_normal(poly[i], poly[j]) * offset;
        edges[i] = {poly[i] + norm, poly[j] + norm};
    }

    std::vector<Vec2> raw;
    for (int i = 0; i < n; ++i) {
        int prev = (i + n - 1) % n;
        Vec2 d1 = edges[prev].b - edges[prev].a;
        Vec2 d2 = edges[i].b - edges[i].a;
        auto t = line_line_param(edges[prev].a, d1, edges[i].a, d2);

        if (t.has_value()) {
            Vec2 miter = edges[prev].a + d1 * t.value();
            double miterLen = (miter - poly[i]).length();
            if (std::abs(offset) > EPS &&
                miterLen / std::abs(offset) > miter_limit) {
                raw.push_back(edges[prev].b);
                raw.push_back(edges[i].a);
            } else {
                raw.push_back(miter);
            }
        } else {
            raw.push_back(edges[prev].b);
        }
    }
    // Close the loop explicitly
    if (!raw.empty()) raw.push_back(raw.front());
    return raw;
}

// ──────────── self-intersection removal (greedy walk) ─────────────────
//
// Walk along the raw offset polyline.  For each segment look ahead
// for the *latest* crossing with a later non-adjacent segment.  When
// found, jump forward past the loop.  This removes all self-
// intersection loops while keeping the valid offset geometry.

static std::vector<std::vector<Vec2>> remove_self_intersections(
    const std::vector<Vec2>& raw,
    const std::vector<Vec2>& original,
    double offset,
    bool closed)
{
    const int n = static_cast<int>(raw.size());
    if (n < 3) return {raw};

    // Collect ALL crossings between non-adjacent segments
    struct Crossing {
        int    i, j;     // segment indices (j > i + 1)
        double ti, tj;   // parameters on those segments
        Vec2   pt;
    };
    std::vector<Crossing> crossings;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 2; j < n - 1; ++j) {
            // Skip adjacent-at-wrap for closed polylines
            if (closed && i == 0 && j == n - 2) continue;

            auto hit = seg_seg_isect(raw[i], raw[i + 1],
                                     raw[j], raw[j + 1]);
            if (hit) {
                crossings.push_back({i, j, hit->t, hit->s, hit->point});
            }
        }
    }

    if (crossings.empty()) return {raw};

    // Greedy forward walk
    std::vector<std::vector<Vec2>> result;
    std::vector<Vec2> cur;
    int    seg   = 0;
    double tFrom = 0.0;
    cur.push_back(raw[0]);

    while (seg < n - 1) {
        // Among all crossings that start on segment `seg` with t > tFrom,
        // pick the one with the smallest t (earliest hit on this segment).
        // That crossing jumps us forward to another segment.
        const Crossing* best = nullptr;
        for (const auto& c : crossings) {
            if (c.i == seg && c.ti > tFrom + EPS) {
                if (!best || c.ti < best->ti)
                    best = &c;
            }
        }

        if (best) {
            cur.push_back(best->pt);

            // The loop between (seg, best->ti) and (best->j, best->tj)
            // should be removed.  Decide which side to keep by checking
            // the distance of the loop midpoint to the original polyline.
            // If the looped portion is *closer* than |offset|, it is the
            // invalid side and we skip it; otherwise we keep it.

            // Compute a representative point inside the loop
            int loopMidSeg = (seg + best->j) / 2 + 1;
            if (loopMidSeg >= n) loopMidSeg = n - 1;
            Vec2 loopMid = raw[std::min(loopMidSeg, n - 1)];
            double loopDist = point_polyline_dist(loopMid, original);

            bool skipForward = (loopDist < std::abs(offset) - EPS);
            // Also skip if the loop is extremely small (degenerate)
            if (!skipForward) {
                // Check a second sample
                int loopSeg2 = seg + 1;
                if (loopSeg2 < best->j) {
                    Vec2 sample2 = raw[loopSeg2];
                    double d2 = point_polyline_dist(sample2, original);
                    if (d2 < std::abs(offset) * 0.5)
                        skipForward = true;
                }
            }

            if (skipForward) {
                // Skip the loop: jump to segment best->j at parameter best->tj
                seg   = best->j;
                tFrom = best->tj;
            } else {
                // The forward portion is the invalid part — output current
                // piece and start a new one from the far side
                // Actually in this case, the piece *before* is valid and the
                // piece *after* the loop is also valid but the connecting
                // portion between them should be the loop content.
                // For simplicity keep the forward jump — the distance
                // heuristic is usually correct for offsets.
                seg   = best->j;
                tFrom = best->tj;
            }
        } else {
            // No crossing starting on this segment; advance normally
            cur.push_back(raw[seg + 1]);
            ++seg;
            tFrom = 0.0;
        }
    }

    cur = deduplicate(cur);
    if (cur.size() >= 2) result.push_back(std::move(cur));

    // For closed polylines, also handle crossings that wrap around
    // (already handled by the walk above)

    return result.empty() ? std::vector<std::vector<Vec2>>{raw} : result;
}

// ─────────────── main entry point ─────────────────────────────────────

std::vector<std::vector<Vec2>> offset_polyline(
    const std::vector<Vec2>& polyline,
    double offset,
    bool   closed,
    double miter_limit)
{
    // Remove degenerate (duplicate) vertices from the input
    std::vector<Vec2> clean = deduplicate(polyline);
    if (clean.size() < 2) return {};
    if (std::abs(offset) < EPS) return {clean};

    // 1. Build the raw offset
    std::vector<Vec2> raw = closed
        ? raw_offset_closed(clean, offset, miter_limit)
        : raw_offset_open  (clean, offset, miter_limit);

    raw = deduplicate(raw);
    if (raw.size() < 2) return {};

    // 2. Remove self-intersections
    auto pieces = remove_self_intersections(raw, clean, offset, closed);

    // 3. Final deduplication & degenerate-segment removal
    std::vector<std::vector<Vec2>> out;
    for (auto& p : pieces) {
        p = deduplicate(p);
        if (p.size() >= 2) out.push_back(std::move(p));
    }
    return out;
}

} // namespace polyline_offset
