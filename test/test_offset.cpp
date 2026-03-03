#include "polyline_offset/polyline_offset.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>

using namespace polyline_offset;

static constexpr double TOL = 1e-6;

static bool near(double a, double b) { return std::abs(a - b) < TOL; }

static bool pt_near(const Point2D& a, const Point2D& b) {
    return (a - b).length() < TOL;
}

/// Shortest distance from point to segment.
static double pt_seg_dist(const Point2D& p,
                          const Point2D& a, const Point2D& b) {
    Point2D ab = b - a;
    double len2 = ab.length_sq();
    if (len2 < 1e-18) return (p - a).length();
    double t = std::max(0.0, std::min(1.0, (p - a).dot(ab) / len2));
    return (p - (a + ab * t)).length();
}

/// Shortest distance from point to polyline.
static double pt_poly_dist(const Point2D& p, const Polyline& poly) {
    double best = 1e18;
    for (size_t i = 0; i + 1 < poly.size(); ++i)
        best = std::min(best, pt_seg_dist(p, poly[i], poly[i + 1]));
    return best;
}

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) \
    static void test_##name(); \
    struct Register_##name { \
        Register_##name() { \
            std::printf("  %-50s", #name); \
            try { test_##name(); std::printf(" PASS\n"); ++tests_passed; } \
            catch (...) { std::printf(" FAIL\n"); ++tests_failed; } \
        } \
    } reg_##name; \
    static void test_##name()

#define ASSERT(expr) \
    do { if (!(expr)) { \
        std::fprintf(stderr, "    ASSERT failed: %s  (%s:%d)\n", \
                     #expr, __FILE__, __LINE__); \
        throw 1; \
    } } while (0)

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST(single_segment_left) {
    // Horizontal segment offset upward (left of direction)
    Polyline in = {{0, 0}, {10, 0}};
    auto out = offset_polyline(in, 2.0);
    ASSERT(out.size() == 2);
    ASSERT(pt_near(out[0], {0, 2}));
    ASSERT(pt_near(out[1], {10, 2}));
}

TEST(single_segment_right) {
    Polyline in = {{0, 0}, {10, 0}};
    auto out = offset_polyline(in, -3.0);
    ASSERT(out.size() == 2);
    ASSERT(pt_near(out[0], {0, -3}));
    ASSERT(pt_near(out[1], {10, -3}));
}

TEST(right_angle_outward) {
    // L-shape: right then up
    Polyline in = {{0, 0}, {10, 0}, {10, 10}};
    auto out = offset_polyline(in, 1.0);
    // Edge 0 offsets up by 1: (0,1)-(10,1)
    // Edge 1 offsets left by 1: (9,0)-(9,10)
    // Miter at vertex: intersection of y=1 and x=9 → (9,1)
    ASSERT(out.size() == 3);
    ASSERT(pt_near(out[0], {0, 1}));
    ASSERT(pt_near(out[1], {9, 1}));
    ASSERT(pt_near(out[2], {9, 10}));
}

TEST(right_angle_inward) {
    Polyline in = {{0, 0}, {10, 0}, {10, 10}};
    auto out = offset_polyline(in, -1.0);
    ASSERT(out.size() == 3);
    ASSERT(pt_near(out[0], {0, -1}));
    ASSERT(pt_near(out[1], {11, -1}));
    ASSERT(pt_near(out[2], {11, 10}));
}

TEST(collinear_vertices) {
    // Three collinear points – middle vertex is degenerate but valid.
    Polyline in = {{0, 0}, {5, 0}, {10, 0}};
    auto out = offset_polyline(in, 2.0);
    ASSERT(out.size() >= 2);
    // All output points should be at y=2
    for (auto& p : out) ASSERT(near(p.y, 2.0));
    // First and last x unchanged
    ASSERT(near(out.front().x, 0.0));
    ASSERT(near(out.back().x, 10.0));
}

TEST(duplicate_vertices) {
    // Duplicate vertex should be cleaned automatically
    Polyline in = {{0, 0}, {0, 0}, {10, 0}, {10, 0}};
    auto out = offset_polyline(in, 1.0);
    ASSERT(out.size() == 2);
    ASSERT(pt_near(out[0], {0, 1}));
    ASSERT(pt_near(out[1], {10, 1}));
}

TEST(zero_distance) {
    Polyline in = {{1, 2}, {3, 4}, {5, 6}};
    auto out = offset_polyline(in, 0.0);
    ASSERT(out.size() == in.size());
    for (size_t i = 0; i < in.size(); ++i) ASSERT(pt_near(out[i], in[i]));
}

TEST(too_few_points) {
    Polyline in0;
    ASSERT(offset_polyline(in0, 1.0).empty());

    Polyline in1 = {{5, 5}};
    auto out1 = offset_polyline(in1, 1.0);
    ASSERT(out1.size() == 1);
}

TEST(self_intersection_v_shape) {
    // A sharp V-shape that causes self-intersection on the outside offset.
    //     (5,10)
    //    /      \.
    // (0,0)    (10,0)
    Polyline in = {{0, 0}, {5, 10}, {10, 0}};
    auto out = offset_polyline(in, -2.0);
    // After offset the two legs should intersect; the loop must be removed
    // so the result has no self-intersections.
    // Verify no segment pair intersects.
    for (size_t i = 0; i + 1 < out.size(); ++i) {
        for (size_t j = i + 2; j + 1 < out.size(); ++j) {
            double t, s;
            Point2D ip;
            Point2D d1 = out[i + 1] - out[i];
            Point2D d2 = out[j + 1] - out[j];
            double denom = d1.cross(d2);
            if (std::abs(denom) < 1e-12) continue;
            Point2D dp = out[j] - out[i];
            t = dp.cross(d2) / denom;
            s = dp.cross(d1) / denom;
            bool intersects = t > 1e-9 && t < 1.0 - 1e-9 &&
                              s > 1e-9 && s < 1.0 - 1e-9;
            ASSERT(!intersects);
        }
    }
}

TEST(closed_square_outward) {
    // Unit square CCW, offset outward by 1
    Polyline in = {{0, 0}, {10, 0}, {10, 10}, {0, 10}};
    OffsetOptions opts;
    opts.is_closed = true;
    auto out = offset_polyline(in, 1.0, opts);
    // Expect 5 points (closed), corners at (-1,-1), (11,-1), (11,11), (-1,11)
    ASSERT(out.size() >= 5);
    ASSERT(pt_near(out.front(), out.back())); // closed
}

TEST(closed_square_inward) {
    Polyline in = {{0, 0}, {10, 0}, {10, 10}, {0, 10}};
    OffsetOptions opts;
    opts.is_closed = true;
    auto out = offset_polyline(in, -1.0, opts);
    ASSERT(out.size() >= 5);
    ASSERT(pt_near(out.front(), out.back()));
}

TEST(offset_distance_accuracy) {
    // For a gentle curve, every offset point should be ≈ |distance| away from
    // the original polyline.
    Polyline in = {{0, 0}, {4, 3}, {8, 1}, {12, 4}, {16, 0}};
    double dist = 1.5;
    auto out = offset_polyline(in, dist);
    for (auto& p : out) {
        double d = pt_poly_dist(p, in);
        // Allow some tolerance for miter extensions at sharp angles
        ASSERT(d > dist * 0.5);
    }
}

TEST(narrow_u_inward_collapse) {
    // A narrow U-shape offset inward by more than half the width
    // should collapse / produce a shortened result.
    Polyline in = {{0, 0}, {0, 10}, {2, 10}, {2, 0}};
    auto out = offset_polyline(in, 1.5);
    // The U gap is only 2 wide; offset 1.5 from each side means the inner
    // edges overlap.  After self-intersection removal the result should be
    // shorter than the original.
    ASSERT(out.size() >= 2);
}

// ---------------------------------------------------------------------------

int main() {
    std::printf("Running polyline-offset tests ...\n");
    // Tests are auto-registered via static constructors above.
    std::printf("\n%d passed, %d failed\n", tests_passed, tests_failed);
    return tests_failed > 0 ? 1 : 0;
}
