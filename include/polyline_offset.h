#pragma once

#include <cmath>
#include <vector>
#include <optional>
#include <algorithm>

namespace polyline_offset {

/// 浮点数比较所用的容差（精度阈值）。
constexpr double EPS = 1e-9;

// ───────────────────────────── 2-D vector ─────────────────────────────
// 二维向量结构体，提供常用的向量运算。
struct Vec2 {
    double x = 0, y = 0;

    Vec2() = default;
    Vec2(double x, double y) : x(x), y(y) {}

    Vec2 operator+(const Vec2& o) const { return {x + o.x, y + o.y}; }  // 向量加法
    Vec2 operator-(const Vec2& o) const { return {x - o.x, y - o.y}; }  // 向量减法
    Vec2 operator*(double s)      const { return {x * s, y * s}; }       // 标量乘法
    // 近似相等：两分量之差均在 EPS 以内视为相等
    bool operator==(const Vec2& o) const {
        return std::abs(x - o.x) < EPS && std::abs(y - o.y) < EPS;
    }

    double dot  (const Vec2& o) const { return x * o.x + y * o.y; }     // 点积
    double cross(const Vec2& o) const { return x * o.y - y * o.x; }     // 叉积（2D 标量）
    double length()   const { return std::sqrt(x * x + y * y); }        // 向量模长
    double lengthSq() const { return x * x + y * y; }                   // 模长的平方（避免开方）
    Vec2   perp()     const { return {-y, x}; }          // 左法向量（逆时针旋转 90°）

    // 返回单位向量；若模长过小则返回零向量
    Vec2 normalized() const {
        double l = length();
        return l > 1e-12 ? Vec2{x / l, y / l} : Vec2{0, 0};
    }
};

// ──────────────────────── segment intersection ────────────────────────
// 线段交点信息结构体。
struct SegIsect {
    double t;        // 第一条线段上的参数 t，取值范围 [0,1]
    double s;        // 第二条线段上的参数 s，取值范围 [0,1]
    Vec2   point;    // 交点坐标
};

/// 直线–直线求交（返回第一条直线上的参数 t）。
/// P + t*D 与 Q + s*E 的交点参数；若平行则返回 nullopt。
std::optional<double> line_line_param(
    const Vec2& P, const Vec2& D,
    const Vec2& Q, const Vec2& E);

/// 线段–线段求交；不相交时返回 nullopt。
std::optional<SegIsect> seg_seg_isect(
    const Vec2& A, const Vec2& B,
    const Vec2& C, const Vec2& D);

/// 计算点到折线的最短距离。
double point_polyline_dist(
    const Vec2& P,
    const std::vector<Vec2>& poly);

// ──────────────────── main offset function ────────────────────────────
/// 对开放或封闭折线进行偏移操作。
///
/// @param polyline   折线的有序顶点列表。
/// @param offset     偏移距离（正值 = 向左/向外偏移，
///                   负值 = 向右/向内偏移）。
/// @param closed     折线是否封闭（即首尾相连构成环）。
/// @param miter_limit 斜切长度与偏移距离之比的上限。
///                    超过此值的顶点将改用斜角（bevel）连接。
/// @return 一条或多条偏移后的折线（自相交消除可能导致结果被分割成多段）。
std::vector<std::vector<Vec2>> offset_polyline(
    const std::vector<Vec2>& polyline,
    double offset,
    bool   closed      = false,
    double miter_limit  = 2.0);

} // namespace polyline_offset
