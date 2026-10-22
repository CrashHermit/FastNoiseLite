## Getting Started

Here's an example for creating a 128x128 array of OpenSimplex2 noise

```cpp
// Create and configure FastNoise object
FastNoiseLite noise;
noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

// Gather noise data
std::vector<float> noiseData(128 * 128);
int index = 0;

for (int y = 0; y < 128; y++)
{
    for (int x = 0; x < 128; x++)
    {
        noiseData[index++] = noise.GetNoise((float)x, (float)y);
    }
}

// Do something with this data...
```

## Analytic OpenSimplex2S gradients (C++ fork extension)

This fork adds `FastNoiseLite::NoiseGradient3D`, containing `value`, `dx`, `dy`,
and `dz`, and `GetOpenSimplex2SWithGradient(x, y, z)`. Both float and double
coordinates are accepted; the kernel and returned components use float
arithmetic.

```cpp
FastNoiseLite noise(123);
noise.SetFrequency(1.0f);
auto sample = noise.GetOpenSimplex2SWithGradient(0.2, -0.4, 0.8);
// sample.value: noise; sample.dx/dy/dz: derivatives in the input coordinates.
```

The new method always samples raw 3D OpenSimplex2S. It honors seed, frequency,
and `RotationType3D`, independently of the selected noise type or fractal mode.
Its value matches `GetNoise` configured with `NoiseType_OpenSimplex2S` and
`FractalType_None`. Internal rotation and frequency are included in the returned
derivatives. Floating-point coordinate quantization is not differentiated.

Value and gradient share the corrected lattice traversal, hashes, gradient
tables, and normalization. Four opposite-corner pairs on each of the two
cubic lattices supply eight candidates in total. This covers every positive
contribution: opposite vertices are twice the support radius apart, so the
unselected member cannot have positive support. At a pair-selection boundary
both contributions are zero.

Value-only calls instantiate the same traversal without gradient calculations
or gradient allocations. No finite differences are used in production. This
API provides first derivatives only; it does not add gradient APIs to fractals,
domain warps, other algorithms, or other language ports.

Run the standalone C++ regression check from the repository root:

```sh
c++ -std=c++11 -O2 -fwrapv Cpp/tests/opensimplex2s_gradients.cpp -o /tmp/fnl-gradients
/tmp/fnl-gradients
```

The check compares fused and value-only outputs exactly, checks derivatives
numerically across seeds, orientations and frequencies, verifies the former
seams, and compares against an independent full sixteen-vertex support sum.
Numerical differentiation is used only by the check.

The previous three-dimensional OpenSimplex2S implementation omitted some
positive corner contributions through conditional shortcuts. This correction
intentionally changes scalar values where those contributions were missing;
floating-point accumulation order can also change. The existing seed hashing,
gradient table, half-lattice index convention and normalization are retained.
Two-dimensional OpenSimplex2S and the other noise algorithms are unchanged.

The underlying compact-support polynomial contributions and their first
derivatives join continuously in real arithmetic. Finite-precision coordinate
quantization and rounding remain; this is not a promise of unlimited numerical
precision or of differentiability of downstream consumers.
