# polyline-offset

高效多段线 Offset 算法（C++）——支持外扩/内缩、顶点退化、自相交处理

## Features

- **Outward / Inward offset** – positive offset shifts left (outward for CCW polygons), negative shifts right (inward).
- **Vertex degeneration** – duplicate and near-duplicate vertices are automatically removed.
- **Self-intersection removal** – loops caused by large offsets on concave geometry are detected and trimmed.
- **Miter-limit bevel** – sharp corners exceeding the configurable miter limit are beveled.
- **Open & closed polylines** – works with both open paths and closed polygons.

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
```

## Run tests

```bash
cd build
ctest --output-on-failure
```

## API

```cpp
#include "polyline_offset.h"
using namespace polyline_offset;

// Offset an open polyline to the left by 1.0
std::vector<Vec2> poly = {{0,0}, {10,0}, {10,10}};
auto result = offset_polyline(poly, 1.0);

// Offset a closed polygon inward by 0.5
std::vector<Vec2> square = {{0,0}, {10,0}, {10,10}, {0,10}};
auto result2 = offset_polyline(square, -0.5, /*closed=*/true);

// Custom miter limit (default is 2.0)
auto result3 = offset_polyline(poly, 1.0, false, 4.0);
```

The function returns `std::vector<std::vector<Vec2>>` – one or more polylines, because self-intersection removal can split the result into separate pieces.