#include "polyline_offset.h"
#include <cassert>
#include <cmath>
#include <limits>

namespace polyline_offset {

// ──────────────────────── helper geometry ─────────────────────────────
// 辅助几何计算函数

// 求两条直线的交点参数 t（第一条直线 P + t*D）。
// 若两直线平行（叉积接近零）则返回 nullopt。
std::optional<double> line_line_param(
    const Vec2& P, const Vec2& D,
    const Vec2& Q, const Vec2& E)
{
    double denom = D.cross(E);
    if (std::abs(denom) < EPS) return std::nullopt;  // 平行，无交点
    return (Q - P).cross(E) / denom;
}

// 求两条线段的交点。
// 通过参数方程 A+t*AB 和 C+s*CD 求解 t 和 s；
// 若 t、s 均在 [0,1] 范围内则存在交点，否则返回 nullopt。
std::optional<SegIsect> seg_seg_isect(
    const Vec2& A, const Vec2& B,
    const Vec2& C, const Vec2& D)
{
    Vec2   AB = B - A, CD = D - C;
    double denom = AB.cross(CD);
    if (std::abs(denom) < EPS) return std::nullopt;  // 平行或共线，无交点

    Vec2   AC = C - A;
    double t  = AC.cross(CD) / denom;
    double s  = AC.cross(AB) / denom;

    // t、s 超出 [0,1] 范围说明交点不在线段上
    if (t < -EPS || t > 1 + EPS || s < -EPS || s > 1 + EPS)
        return std::nullopt;

    t = std::clamp(t, 0.0, 1.0);
    s = std::clamp(s, 0.0, 1.0);
    return SegIsect{t, s, A + AB * t};
}

// 计算点 P 到折线的最短距离。
// 遍历折线的每条线段，对每段求点到线段的距离，取最小值。
double point_polyline_dist(const Vec2& P, const std::vector<Vec2>& poly)
{
    double best = std::numeric_limits<double>::max();
    for (size_t i = 0; i + 1 < poly.size(); ++i) {
        Vec2   AB   = poly[i + 1] - poly[i];
        double len2 = AB.lengthSq();
        // 将 P 投影到线段 AB 上，参数 t 夹紧到 [0,1] 以限制在线段范围内
        double t    = (len2 > 1e-24)
                        ? std::clamp((P - poly[i]).dot(AB) / len2, 0.0, 1.0)
                        : 0.0;
        double d = (P - (poly[i] + AB * t)).length();
        best = std::min(best, d);
    }
    return best;
}

// ──────────── remove near-duplicate consecutive vertices ──────────────
// 去除连续的近重复顶点（顶点去重）。
// 当相邻两点距离的平方小于 tol² 时视为重复并跳过。

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
// 计算线段 AB 的单位左法向量（即垂直于 AB 且朝左的单位向量）。

static Vec2 seg_normal(const Vec2& A, const Vec2& B)
{
    return (B - A).perp().normalized();
}

// ─────────────── build raw offset polyline (open case) ────────────────
//
// 构建开放折线的原始偏移折线。
//
// 对每条原始边沿其左法向量方向移动 `offset` 距离，得到偏移边。
// 相邻偏移边通过求交点（斜切连接）进行连接；
// 若斜切比超过 miter_limit，则改用斜角（bevel）连接——插入当前边的终点
// 和下一条边的起点，形成一段短的斜切线段。

static std::vector<Vec2> raw_offset_open(
    const std::vector<Vec2>& poly,
    double offset,
    double miter_limit)
{
    const int n = static_cast<int>(poly.size());
    if (n < 2) return poly;

    // 预计算各偏移边：第 i 条边偏移后为 (edges[i].a, edges[i].b)
    struct Edge { Vec2 a, b; };
    std::vector<Edge> edges(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        Vec2 norm = seg_normal(poly[i], poly[i + 1]) * offset;
        edges[i] = {poly[i] + norm, poly[i + 1] + norm};
    }

    std::vector<Vec2> raw;
    raw.push_back(edges[0].a);  // 起点直接取第一条偏移边的起点

    for (int i = 0; i < n - 2; ++i) {
        // 求第 i 条偏移边与第 i+1 条偏移边的直线交点
        Vec2 d1 = edges[i].b - edges[i].a;
        Vec2 d2 = edges[i + 1].b - edges[i + 1].a;
        auto t = line_line_param(edges[i].a, d1, edges[i + 1].a, d2);

        if (t.has_value()) {
            Vec2 miter = edges[i].a + d1 * t.value();
            double miterLen = (miter - poly[i + 1]).length();
            if (std::abs(offset) > EPS &&
                miterLen / std::abs(offset) > miter_limit) {
                // 斜切比超限，改用斜角连接：插入当前边终点和下一边起点
                raw.push_back(edges[i].b);
                raw.push_back(edges[i + 1].a);
            } else {
                // 斜切比在限制范围内，使用斜切交点
                raw.push_back(miter);
            }
        } else {
            // 相邻偏移边平行，直接使用当前边的终点
            raw.push_back(edges[i].b);
        }
    }

    raw.push_back(edges.back().b);  // 终点取最后一条偏移边的终点
    return raw;
}

// ─────────────── build raw offset polyline (closed case) ──────────────
// 构建封闭折线（多边形）的原始偏移折线。
// 与开放情形类似，但需处理首尾边之间的连接，并在末尾显式闭合结果。

static std::vector<Vec2> raw_offset_closed(
    const std::vector<Vec2>& poly,
    double offset,
    double miter_limit)
{
    int n = static_cast<int>(poly.size());
    // 若输入已显式闭合（首尾顶点相同），去掉重复的末尾顶点
    if (n >= 2 && poly.front() == poly.back()) --n;
    if (n < 2) return poly;

    struct Edge { Vec2 a, b; };
    std::vector<Edge> edges(n);
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;  // 下一顶点索引（循环）
        Vec2 norm = seg_normal(poly[i], poly[j]) * offset;
        edges[i] = {poly[i] + norm, poly[j] + norm};
    }

    std::vector<Vec2> raw;
    for (int i = 0; i < n; ++i) {
        int prev = (i + n - 1) % n;  // 前一条边的索引（循环）
        Vec2 d1 = edges[prev].b - edges[prev].a;
        Vec2 d2 = edges[i].b - edges[i].a;
        auto t = line_line_param(edges[prev].a, d1, edges[i].a, d2);

        if (t.has_value()) {
            Vec2 miter = edges[prev].a + d1 * t.value();
            double miterLen = (miter - poly[i]).length();
            if (std::abs(offset) > EPS &&
                miterLen / std::abs(offset) > miter_limit) {
                // 斜切比超限，改用斜角连接
                raw.push_back(edges[prev].b);
                raw.push_back(edges[i].a);
            } else {
                raw.push_back(miter);
            }
        } else {
            // 相邻偏移边平行
            raw.push_back(edges[prev].b);
        }
    }
    // 显式闭合环：将首顶点追加到末尾
    if (!raw.empty()) raw.push_back(raw.front());
    return raw;
}

// ──────────── self-intersection removal (greedy walk) ─────────────────
//
// 自相交消除（贪心前向遍历）。
//
// 沿原始偏移折线向前遍历，对每条线段向后查找与非相邻线段的最早交叉点。
// 发现交叉后，跳过中间形成的自相交环，保留有效的偏移几何形状。
// 判断哪一侧为无效环的依据是：检查环内采样点到原始折线的距离；
// 若该距离小于偏移量绝对值，说明该部分位于偏移带内侧，应被跳过。

static std::vector<std::vector<Vec2>> remove_self_intersections(
    const std::vector<Vec2>& raw,
    const std::vector<Vec2>& original,
    double offset,
    bool closed)
{
    const int n = static_cast<int>(raw.size());
    if (n < 3) return {raw};

    // 收集所有非相邻线段之间的交叉点
    struct Crossing {
        int    i, j;     // 线段索引（j > i + 1）
        double ti, tj;   // 各自线段上的参数值
        Vec2   pt;       // 交叉点坐标
    };
    std::vector<Crossing> crossings;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 2; j < n - 1; ++j) {
            // 对封闭折线跳过首尾相邻的边对（避免误判）
            if (closed && i == 0 && j == n - 2) continue;

            auto hit = seg_seg_isect(raw[i], raw[i + 1],
                                     raw[j], raw[j + 1]);
            if (hit) {
                crossings.push_back({i, j, hit->t, hit->s, hit->point});
            }
        }
    }

    if (crossings.empty()) return {raw};  // 无自相交，直接返回

    // 贪心前向遍历：逐段走，遇到交叉点则跳过自相交环
    std::vector<std::vector<Vec2>> result;
    std::vector<Vec2> cur;
    int    seg   = 0;
    double tFrom = 0.0;
    cur.push_back(raw[0]);

    while (seg < n - 1) {
        // 在当前线段 seg 上，找参数值大于 tFrom 的最早交叉点
        const Crossing* best = nullptr;
        for (const auto& c : crossings) {
            if (c.i == seg && c.ti > tFrom + EPS) {
                if (!best || c.ti < best->ti)
                    best = &c;
            }
        }

        if (best) {
            cur.push_back(best->pt);  // 将交叉点加入当前折线段

            // 判断是否应跳过从当前交叉点到 best->j 的自相交环。
            // 计算环内代表点（取环中间线段上的点）到原始折线的距离：
            // 若距离小于偏移量绝对值，说明该环位于无效区域，应跳过。

            // 计算环内中间线段的索引
            int loopMidSeg = seg + (best->j - seg) / 2 + 1;
            if (loopMidSeg >= n - 1) loopMidSeg = n - 2;
            if (loopMidSeg < 0)      loopMidSeg = 0;
            Vec2 loopMid = raw[loopMidSeg];
            double loopDist = point_polyline_dist(loopMid, original);

            bool skipForward = (loopDist < std::abs(offset) - EPS);
            // 补充检查：若第二个采样点距原始折线不足偏移量的一半，
            // 也视为无效环（启发式阈值 0.5）。
            constexpr double kLoopDistThreshold = 0.5;
            if (!skipForward) {
                // 检查环起始处的第二个采样点
                int loopSeg2 = seg + 1;
                if (loopSeg2 < n && loopSeg2 < best->j) {
                    Vec2 sample2 = raw[loopSeg2];
                    double d2 = point_polyline_dist(sample2, original);
                    if (d2 < std::abs(offset) * kLoopDistThreshold)
                        skipForward = true;
                }
            }

            if (skipForward) {
                // 跳过自相交环：直接跳到 best->j 线段的 best->tj 处继续
                seg   = best->j;
                tFrom = best->tj;
            } else {
                // 启发式判断前向部分为无效段：同样跳过
                // （距离启发式通常对偏移情形正确）
                seg   = best->j;
                tFrom = best->tj;
            }
        } else {
            // 当前线段无交叉点，正常前进到下一线段
            cur.push_back(raw[seg + 1]);
            ++seg;
            tFrom = 0.0;
        }
    }

    cur = deduplicate(cur);
    if (cur.size() >= 2) result.push_back(std::move(cur));

    // 封闭折线的环绕交叉已在上述遍历中处理

    return result.empty() ? std::vector<std::vector<Vec2>>{raw} : result;
}

// ─────────────── main entry point ─────────────────────────────────────
// 主入口函数：对折线执行偏移操作，返回一条或多条偏移后的折线。

std::vector<std::vector<Vec2>> offset_polyline(
    const std::vector<Vec2>& polyline,
    double offset,
    bool   closed,
    double miter_limit)
{
    // 第 0 步：去除输入中的退化（重复）顶点
    std::vector<Vec2> clean = deduplicate(polyline);
    if (clean.size() < 2) return {};      // 少于两点无法构成折线
    if (std::abs(offset) < EPS) return {clean};  // 偏移量为零，直接返回清理后的输入

    // 第 1 步：构建原始偏移折线
    std::vector<Vec2> raw = closed
        ? raw_offset_closed(clean, offset, miter_limit)
        : raw_offset_open  (clean, offset, miter_limit);

    raw = deduplicate(raw);
    if (raw.size() < 2) return {};

    // 第 2 步：消除自相交
    auto pieces = remove_self_intersections(raw, clean, offset, closed);

    // 第 3 步：对结果进行最终去重，并过滤掉退化（点数不足）的段
    std::vector<std::vector<Vec2>> out;
    for (auto& p : pieces) {
        p = deduplicate(p);
        if (p.size() >= 2) out.push_back(std::move(p));
    }
    return out;
}

} // namespace polyline_offset
