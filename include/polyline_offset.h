#pragma once

#include <cmath>
#include <vector>
#include <optional>
#include <algorithm>

namespace polyline_offset {

// ───────────────────────────── 2-D vector ─────────────────────────────
struct Vec2 {
    double x = 0, y = 0;

    Vec2() = default;
    Vec2(double x, double y) : x(x), y(y) {}

    Vec2 operator+(const Vec2& o) const { return {x + o.x, y + o.y}; }
    Vec2 operator-(const Vec2& o) const { return {x - o.x, y - o.y}; }
    Vec2 operator*(double s)      const { return {x * s, y * s}; }
    bool operator==(const Vec2& o) const {
        return std::abs(x - o.x) < 1e-9 && std::abs(y - o.y) < 1e-9;
    }

    double dot  (const Vec2& o) const { return x * o.x + y * o.y; }
    double cross(const Vec2& o) const { return x * o.y - y * o.x; }
    double length()   const { return std::sqrt(x * x + y * y); }
    double lengthSq() const { return x * x + y * y; }
    Vec2   perp()     const { return {-y, x}; }          // left normal

    Vec2 normalized() const {
        double l = length();
        return l > 1e-12 ? Vec2{x / l, y / l} : Vec2{0, 0};
    }
};

// ──────────────────────── segment intersection ────────────────────────
struct SegIsect {
    double t;        // parameter on first segment  [0,1]
    double s;        // parameter on second segment [0,1]
    Vec2   point;
};

/// Line–line intersection (returns parameter t on first line).
std::optional<double> line_line_param(
    const Vec2& P, const Vec2& D,
    const Vec2& Q, const Vec2& E);

/// Segment–segment intersection.
std::optional<SegIsect> seg_seg_isect(
    const Vec2& A, const Vec2& B,
    const Vec2& C, const Vec2& D);

/// Minimum distance from a point to a polyline.
double point_polyline_dist(
    const Vec2& P,
    const std::vector<Vec2>& poly);

// ──────────────────── main offset function ────────────────────────────
/// Offset an open or closed polyline.
///
/// @param polyline   Ordered vertices of the polyline.
/// @param offset     Offset distance (positive = left / outward,
///                   negative = right / inward).
/// @param closed     Whether the polyline is closed (forms a loop).
/// @param miter_limit Maximum ratio of miter length to offset distance.
///                    Vertices where the ratio exceeds this are beveled.
/// @return One or more offset polylines (may split due to
///         self-intersection removal).
std::vector<std::vector<Vec2>> offset_polyline(
    const std::vector<Vec2>& polyline,
    double offset,
    bool   closed      = false,
    double miter_limit  = 2.0);

} // namespace polyline_offset
