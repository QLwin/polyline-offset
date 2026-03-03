#pragma once

#include "types.h"

namespace polyline_offset {

/// Options for polyline offset
struct OffsetOptions {
    /// Whether the polyline is closed (forms a polygon).
    /// If true, an edge from the last vertex back to the first is included.
    bool is_closed = false;
};

/// Offset a polyline by the given distance.
///
/// @param polyline  Input polyline (sequence of 2D points, at least 2).
/// @param distance  Offset distance.
///                  Positive = left-side offset (outward for CCW winding).
///                  Negative = right-side offset (inward for CCW winding).
/// @param options   Additional options (e.g. closed polyline).
/// @return The offset polyline with self-intersections removed and degenerate
///         vertices cleaned up.  Returns an empty polyline when the offset
///         causes the geometry to collapse completely.
Polyline offset_polyline(const Polyline& polyline, double distance,
                         const OffsetOptions& options = {});

} // namespace polyline_offset
