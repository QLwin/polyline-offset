#pragma once

#include <cmath>
#include <vector>

namespace polyline_offset {

/// 2D point
struct Point2D {
    double x, y;

    Point2D() : x(0), y(0) {}
    Point2D(double x, double y) : x(x), y(y) {}

    Point2D operator+(const Point2D& p) const { return {x + p.x, y + p.y}; }
    Point2D operator-(const Point2D& p) const { return {x - p.x, y - p.y}; }
    Point2D operator*(double s) const { return {x * s, y * s}; }

    double dot(const Point2D& p) const { return x * p.x + y * p.y; }
    double cross(const Point2D& p) const { return x * p.y - y * p.x; }
    double length() const { return std::sqrt(x * x + y * y); }
    double length_sq() const { return x * x + y * y; }

    Point2D normalized() const {
        double len = length();
        if (len < 1e-12) return {0, 0};
        return {x / len, y / len};
    }
};

/// A polyline is a sequence of 2D points
using Polyline = std::vector<Point2D>;

} // namespace polyline_offset
