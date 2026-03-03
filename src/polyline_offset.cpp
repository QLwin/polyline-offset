#include "polyline_offset/polyline_offset.h"

#include <algorithm>
#include <cmath>

namespace polyline_offset {

namespace {

constexpr double EPS = 1e-9;

// ---- geometry helpers -------------------------------------------------------

/// Unit-length left-hand normal of segment (a, b).
Point2D left_normal(const Point2D& a, const Point2D& b) {
    Point2D d = b - a;
    double len = d.length();
    if (len < EPS) return {0, 0};
    return {-d.y / len, d.x / len};
}

/// Intersect two infinite lines: p1 + t*d1 and p2 + s*d2.
/// Returns false when lines are (nearly) parallel.
bool line_intersect(const Point2D& p1, const Point2D& d1,
                    const Point2D& p2, const Point2D& d2,
                    double& t, double& s) {
    double denom = d1.cross(d2);
    if (std::abs(denom) < EPS) return false;
    Point2D dp = p2 - p1;
    t = dp.cross(d2) / denom;
    s = dp.cross(d1) / denom;
    return true;
}

/// Intersect two finite segments a1-b1 and a2-b2.
/// On success writes parameter values t, s ∈ [0,1] and intersection point ip.
bool segment_intersect(const Point2D& a1, const Point2D& b1,
                       const Point2D& a2, const Point2D& b2,
                       double& t, double& s, Point2D& ip) {
    Point2D d1 = b1 - a1;
    Point2D d2 = b2 - a2;
    if (!line_intersect(a1, d1, a2, d2, t, s)) return false;
    if (t < -EPS || t > 1.0 + EPS || s < -EPS || s > 1.0 + EPS) return false;
    t = std::max(0.0, std::min(1.0, t));
    s = std::max(0.0, std::min(1.0, s));
    ip = a1 + d1 * t;
    return true;
}


// ---- degeneration -----------------------------------------------------------

/// Remove consecutive duplicate (or near-duplicate) vertices.
Polyline remove_degenerate_vertices(const Polyline& poly) {
    if (poly.size() < 2) return poly;
    Polyline out;
    out.push_back(poly[0]);
    for (size_t i = 1; i < poly.size(); ++i) {
        if ((poly[i] - out.back()).length() > EPS) {
            out.push_back(poly[i]);
        }
    }
    return out;
}

// ---- raw offset -------------------------------------------------------------

struct OffsetEdge {
    Point2D start, end;
};

/// Compute offset edges (one per original edge) displaced by `dist`.
std::vector<OffsetEdge> build_offset_edges(const Polyline& poly, double dist) {
    std::vector<OffsetEdge> edges;
    for (size_t i = 0; i + 1 < poly.size(); ++i) {
        Point2D n = left_normal(poly[i], poly[i + 1]) * dist;
        edges.push_back({poly[i] + n, poly[i + 1] + n});
    }
    return edges;
}

/// Build the raw offset polyline by joining adjacent offset edges at their
/// line-line intersection (miter join) or falling back to the edge endpoint
/// when edges are parallel.
Polyline join_offset_edges(const std::vector<OffsetEdge>& edges) {
    if (edges.empty()) return {};
    Polyline out;
    out.push_back(edges.front().start);
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
        Point2D d1 = edges[i].end - edges[i].start;
        Point2D d2 = edges[i + 1].end - edges[i + 1].start;
        double t, s;
        if (line_intersect(edges[i].start, d1, edges[i + 1].start, d2, t, s)) {
            out.push_back(edges[i].start + d1 * t);
        } else {
            out.push_back(edges[i].end);
        }
    }
    out.push_back(edges.back().end);
    return out;
}

/// Build the raw offset polyline for a closed polyline.
/// The result is also closed (first == last).
Polyline join_offset_edges_closed(const std::vector<OffsetEdge>& edges) {
    if (edges.empty()) return {};
    size_t n = edges.size();
    Polyline out;
    for (size_t i = 0; i < n; ++i) {
        size_t next = (i + 1) % n;
        Point2D d1 = edges[i].end - edges[i].start;
        Point2D d2 = edges[next].end - edges[next].start;
        double t, s;
        if (line_intersect(edges[i].start, d1, edges[next].start, d2, t, s)) {
            out.push_back(edges[i].start + d1 * t);
        } else {
            out.push_back(edges[i].end);
        }
    }
    out.push_back(out.front()); // close
    return out;
}

// ---- self-intersection removal ----------------------------------------------

/// Iteratively find and remove the first self-intersection loop until none
/// remain.  Each iteration is O(n²); the total number of iterations is bounded
/// by the number of vertices.
Polyline remove_self_intersections(const Polyline& poly) {
    Polyline cur = poly;

    for (int iter = 0; iter < static_cast<int>(poly.size()) * 2 + 10; ++iter) {
        size_t n = cur.size();
        if (n < 3) break;
        size_t nsegs = n - 1;

        // Detect closed polyline (first ≈ last point)
        bool is_closed_poly =
            n >= 4 && (cur.front() - cur.back()).length() < EPS * 10;

        bool found = false;
        for (size_t i = 0; i < nsegs && !found; ++i) {
            for (size_t j = i + 2; j < nsegs; ++j) {
                double t, s;
                Point2D ip;
                if (!segment_intersect(cur[i], cur[i + 1],
                                       cur[j], cur[j + 1], t, s, ip))
                    continue;

                // For closed polylines, segments 0 and nsegs-1 share the
                // closing vertex – ignore that endpoint touch.
                if (is_closed_poly && i == 0 && j == nsegs - 1 &&
                    std::abs(t) < EPS && std::abs(s - 1.0) < EPS)
                    continue;

                // Decide which part to keep.  We keep the "outside" path
                // (the one whose midpoint is closer to the expected offset
                // distance from the original polyline).
                //
                // Path A = cur[0..i] + ip + cur[j+1..end]   (skip the loop)
                // Path B = cur[0..i] + cur[i+1..j] + ip + cur[j+1..end] (keep)
                // We build path A (skip), check its validity, and use it.

                Polyline path;
                for (size_t k = 0; k <= i; ++k) path.push_back(cur[k]);
                path.push_back(ip);
                for (size_t k = j + 1; k < n; ++k) path.push_back(cur[k]);

                cur = remove_degenerate_vertices(path);
                found = true;
                break;
            }
        }
        if (!found) break;
    }
    return cur;
}

} // anonymous namespace

// ---- public API -------------------------------------------------------------

Polyline offset_polyline(const Polyline& input, double distance,
                         const OffsetOptions& options) {
    if (input.size() < 2) return input;
    if (std::abs(distance) < EPS) return input;

    // 1. Clean degenerate vertices
    Polyline poly = remove_degenerate_vertices(input);
    if (poly.size() < 2) return poly;

    // For closed polylines make sure first == last
    bool closed = options.is_closed;
    if (closed) {
        if ((poly.front() - poly.back()).length() > EPS) {
            poly.push_back(poly.front());
        }
    }

    // 2. Build offset edges
    auto edges = build_offset_edges(poly, distance);
    if (edges.empty()) return {};

    // 3. Join offset edges into a raw offset polyline
    Polyline raw;
    if (closed) {
        // For closed polylines, remove the duplicated closing vertex before
        // joining, because join_offset_edges_closed handles the wrap-around.
        std::vector<OffsetEdge> closed_edges = edges;
        // edges already has n-1 edges where n = poly.size() and poly[0]==poly[n-1]
        raw = join_offset_edges_closed(closed_edges);
    } else {
        raw = join_offset_edges(edges);
    }

    // 4. Remove self-intersections
    Polyline cleaned = remove_self_intersections(raw);

    // 5. Final degenerate-vertex cleanup
    cleaned = remove_degenerate_vertices(cleaned);

    return cleaned;
}

} // namespace polyline_offset
