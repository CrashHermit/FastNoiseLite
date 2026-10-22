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

## Analytic OpenSimplex2S derivatives (C++ fork extension)

The C++ fork adds `NoiseGradient3D` and `NoiseHessian3D`, returned by
`GetOpenSimplex2SWithGradient(x, y, z)` and
`GetOpenSimplex2SWithHessian(x, y, z)`. The Hessian is stored as
`sample.hessian[row][column]`, in input-coordinate order. Float and double
coordinates are accepted; the kernel and returned components use float
arithmetic.

Both APIs sample raw 3D OpenSimplex2S with the configured seed, frequency, and
`RotationType3D`, regardless of the selected noise type or fractal mode. Values
and gradients match their existing APIs exactly. The gradient and Hessian
include the orientation/frequency transforms; the Hessian applies those
transforms on both axes. Existing value and gradient-only calls do not
instantiate Hessian accumulation. No numerical derivatives are used in
production.

The three-dimensional kernel uses eight support candidates: four opposite
corner pairs on each of two offset cubic lattices. This includes every positive
contribution while preserving the established hash/index mapping and
normalization. The compact-support contribution Hessians join continuously
across support boundaries in real arithmetic.

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
