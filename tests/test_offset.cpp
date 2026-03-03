#include "polyline_offset.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>

using namespace polyline_offset;

static int g_pass = 0, g_fail = 0;

static void check(bool cond, const char* msg, int line) {
    if (cond) { ++g_pass; }
    else      { ++g_fail; std::printf("FAIL line %d: %s\n", line, msg); }
}
#define CHECK(c) check((c), #c, __LINE__)

// ───────────────── helper: distance between two points ────────────────
static double dist(const Vec2& a, const Vec2& b) {
    return (a - b).length();
}

// ──────── helper: all points of result are ≈ |offset| from original ───
static bool all_points_at_offset(const std::vector<std::vector<Vec2>>& res,
                                 const std::vector<Vec2>& orig,
                                 double offset,
                                 double tol) {
    double expected = std::abs(offset);
    for (auto& poly : res)
        for (auto& p : poly) {
            double d = point_polyline_dist(p, orig);
            if (std::abs(d - expected) > tol) return false;
        }
    return true;
}

// ══════════════════════════════ tests ══════════════════════════════════

static void test_trivial_input() {
    // Single point → empty result
    CHECK(offset_polyline({{0,0}}, 1.0).empty());
    // Two identical points → empty
    CHECK(offset_polyline({{1,1},{1,1}}, 1.0).empty());
    // Zero offset → returns cleaned input
    auto r = offset_polyline({{0,0},{10,0}}, 0.0);
    CHECK(r.size() == 1);
}

static void test_horizontal_line() {
    //  (0,0)──(10,0)  offset +1 → line at y=1
    auto r = offset_polyline({{0,0},{10,0}}, 1.0);
    CHECK(r.size() == 1);
    CHECK(r[0].size() == 2);
    CHECK(std::abs(r[0][0].y - 1.0) < 1e-6);
    CHECK(std::abs(r[0][1].y - 1.0) < 1e-6);

    // Offset -1 → line at y=-1
    auto r2 = offset_polyline({{0,0},{10,0}}, -1.0);
    CHECK(r2.size() == 1);
    CHECK(std::abs(r2[0][0].y - (-1.0)) < 1e-6);
}

static void test_vertical_line() {
    auto r = offset_polyline({{0,0},{0,10}}, 2.0);
    CHECK(r.size() == 1);
    CHECK(std::abs(r[0][0].x - (-2.0)) < 1e-6);
    CHECK(std::abs(r[0][1].x - (-2.0)) < 1e-6);
}

static void test_right_angle() {
    // L-shape:  (0,0)─(10,0)─(10,10)  offset +1
    // The offset should produce three points with a miter at the corner
    auto r = offset_polyline({{0,0},{10,0},{10,10}}, 1.0);
    CHECK(r.size() == 1);
    CHECK(r[0].size() >= 3);

    // First point should be at y=1 (above horizontal segment)
    CHECK(std::abs(r[0].front().y - 1.0) < 1e-6);

    // Offset -1 (other side)
    auto r2 = offset_polyline({{0,0},{10,0},{10,10}}, -1.0);
    CHECK(r2.size() == 1);
    CHECK(std::abs(r2[0].front().y - (-1.0)) < 1e-6);
}

static void test_acute_angle_bevel() {
    // Very sharp V-shape that should trigger miter-limit bevel
    auto r = offset_polyline({{0,0},{5,0.5},{10,0}}, 1.0, false, 2.0);
    CHECK(r.size() >= 1);
    // With bevel, the join at the apex produces an extra vertex
    CHECK(r[0].size() >= 3);
}

static void test_closed_square() {
    // Unit square offset outward by 1
    std::vector<Vec2> sq = {{0,0},{10,0},{10,10},{0,10}};
    auto r = offset_polyline(sq, 1.0, /*closed=*/true);
    CHECK(r.size() >= 1);
    // All vertices of the result should be ≈1 away from the original
    CHECK(all_points_at_offset(r, sq, 1.0, 0.5));
}

static void test_closed_square_inward() {
    // Inward offset (negative = shrink)
    std::vector<Vec2> sq = {{0,0},{10,0},{10,10},{0,10}};
    auto r = offset_polyline(sq, -1.0, /*closed=*/true);
    CHECK(r.size() >= 1);
    CHECK(all_points_at_offset(r, sq, -1.0, 0.5));
}

static void test_u_shape_self_intersection() {
    // U-shape with a narrow gap – large inward offset triggers
    // self-intersection that must be removed.
    //
    //   (0,10)       (3,10)
    //     │             │
    //     │   (1,2)─(2,2)
    //     │             │
    //   (0,0)       (3,0)
    //
    std::vector<Vec2> u = {{0,0},{0,10},{1,10},{1,2},{2,2},{2,10},{3,10},{3,0}};
    auto r = offset_polyline(u, -0.4);
    // Should produce at least one valid piece after self-intersection removal
    CHECK(!r.empty());
    // No piece should have adjacent duplicate vertices
    for (auto& piece : r) {
        for (size_t i = 0; i + 1 < piece.size(); ++i)
            CHECK(dist(piece[i], piece[i+1]) > 1e-9);
    }
}

static void test_collinear_vertices() {
    // Polyline with redundant collinear vertices
    auto r = offset_polyline({{0,0},{5,0},{10,0}}, 1.0);
    CHECK(r.size() == 1);
    // The two segments are collinear, so offset should still be a single line
    CHECK(std::abs(r[0].front().y - 1.0) < 1e-6);
    CHECK(std::abs(r[0].back().y  - 1.0) < 1e-6);
}

static void test_degenerate_duplicate_vertices() {
    // Polyline with duplicate adjacent vertices (degeneration)
    auto r = offset_polyline({{0,0},{0,0},{10,0},{10,0},{10,0}}, 1.0);
    CHECK(r.size() == 1);
    CHECK(std::abs(r[0].front().y - 1.0) < 1e-6);
}

static void test_closed_triangle() {
    // Equilateral-ish triangle
    std::vector<Vec2> tri = {{0,0},{10,0},{5, 8.66}};
    auto r = offset_polyline(tri, 1.0, true);
    CHECK(r.size() >= 1);
    // Outward offset of a triangle should still be a closed polygon
    // with at least 3 unique vertices (+1 closing)
    CHECK(r[0].size() >= 4);
}

static void test_large_offset_degeneracy() {
    // Offset larger than the feature size → polyline collapses
    // Small triangle with large inward offset
    std::vector<Vec2> tri = {{0,0},{2,0},{1,1}};
    auto r = offset_polyline(tri, -5.0, true);
    // Should return empty or very small result (completely degenerate)
    // Either empty or all points very close together is acceptable
    if (!r.empty()) {
        for (auto& piece : r) {
            // Each piece should be small
            CHECK(piece.size() <= 10);
        }
    }
}

static void test_seg_seg_intersection() {
    // Crossing segments
    auto r = seg_seg_isect({0,0},{10,10},{0,10},{10,0});
    CHECK(r.has_value());
    CHECK(std::abs(r->point.x - 5.0) < 1e-6);
    CHECK(std::abs(r->point.y - 5.0) < 1e-6);

    // Parallel segments → no intersection
    auto r2 = seg_seg_isect({0,0},{10,0},{0,1},{10,1});
    CHECK(!r2.has_value());

    // Non-overlapping segments
    auto r3 = seg_seg_isect({0,0},{1,0},{2,0},{3,0});
    CHECK(!r3.has_value());
}

static void test_point_polyline_dist() {
    std::vector<Vec2> line = {{0,0},{10,0}};
    CHECK(std::abs(point_polyline_dist({5, 3}, line) - 3.0) < 1e-6);
    CHECK(std::abs(point_polyline_dist({0, 0}, line) - 0.0) < 1e-6);
    CHECK(std::abs(point_polyline_dist({-1,0}, line) - 1.0) < 1e-6);
}

static void test_multi_segment_offset() {
    // Zigzag with 5 vertices
    std::vector<Vec2> zigzag = {{0,0},{4,3},{8,0},{12,3},{16,0}};
    auto r = offset_polyline(zigzag, 0.5);
    CHECK(r.size() >= 1);
    CHECK(r[0].size() >= 5);
}

static void test_negative_offset_open() {
    // Simple 3-segment polyline, negative offset
    std::vector<Vec2> poly = {{0,0},{10,0},{10,10},{0,10}};
    auto r = offset_polyline(poly, -1.0);
    CHECK(r.size() >= 1);
    // First point should be at y = -1 (below the first segment)
    CHECK(std::abs(r[0].front().y - (-1.0)) < 1e-6);
}

// ══════════════════════════════ main ══════════════════════════════════

int main() {
    test_trivial_input();
    test_horizontal_line();
    test_vertical_line();
    test_right_angle();
    test_acute_angle_bevel();
    test_closed_square();
    test_closed_square_inward();
    test_u_shape_self_intersection();
    test_collinear_vertices();
    test_degenerate_duplicate_vertices();
    test_closed_triangle();
    test_large_offset_degeneracy();
    test_seg_seg_intersection();
    test_point_polyline_dist();
    test_multi_segment_offset();
    test_negative_offset_open();

    std::printf("\n=== Results: %d passed, %d failed ===\n", g_pass, g_fail);
    return g_fail > 0 ? 1 : 0;
}
