# polyline-offset

高效多段线 Offset 算法（C++）——支持外扩/内缩、顶点退化、自相交处理

## Features

- **Outward / inward offset** – positive distance offsets to the left of the
  polyline direction; negative distance offsets to the right.
- **Vertex degeneration** – duplicate or near-duplicate vertices are
  automatically removed before and after offsetting.
- **Self-intersection removal** – loops created by sharp angles are detected
  and removed automatically.
- **Open & closed polylines** – set `OffsetOptions::is_closed = true` for
  closed polylines (polygons).

## Build

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

## Run tests

```bash
cd build
ctest          # or ./test_offset
```

## Quick example

```cpp
#include "polyline_offset/polyline_offset.h"
#include <cstdio>

int main() {
    using namespace polyline_offset;

    Polyline poly = {{0, 0}, {10, 0}, {10, 10}};
    Polyline result = offset_polyline(poly, 1.0);

    for (auto& p : result)
        std::printf("(%.2f, %.2f)\n", p.x, p.y);

    // Closed polyline (polygon) offset
    Polyline square = {{0, 0}, {10, 0}, {10, 10}, {0, 10}};
    OffsetOptions opts;
    opts.is_closed = true;
    Polyline inward = offset_polyline(square, 1.0, opts);   // shrink (CCW)
    Polyline outward = offset_polyline(square, -1.0, opts);  // expand (CCW)
}
```

## API

```cpp
namespace polyline_offset {

struct Point2D { double x, y; /* operators … */ };
using Polyline = std::vector<Point2D>;

struct OffsetOptions {
    bool is_closed = false;
};

Polyline offset_polyline(const Polyline& polyline, double distance,
                         const OffsetOptions& options = {});
}
```

| Parameter | Description |
|-----------|-------------|
| `polyline` | Input polyline (≥ 2 points). |
| `distance` | Offset distance. Positive = left of walking direction; negative = right. For a CCW polygon with `is_closed = true`, positive shrinks and negative expands. |
| `options.is_closed` | When `true`, an edge from the last vertex back to the first is included and the result is also closed. |