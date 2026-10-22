#include "../FastNoiseLite.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

static void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

template <typename Coordinate>
static float check_point(FastNoiseLite& noise, Coordinate x, Coordinate y, Coordinate z,
    float frequency, double lattice_step = 0.001)
{
    const auto sample = noise.GetOpenSimplex2SWithGradient(x, y, z);
    require(sample.value == noise.GetNoise(x, y, z), "gradient sampler changed scalar value");
    const float gradient[] = {sample.dx, sample.dy, sample.dz};
    if (frequency == 0)
    {
        require(sample.dx == 0 && sample.dy == 0 && sample.dz == 0, "zero frequency must have zero gradient");
        return 0;
    }
    const Coordinate h = (Coordinate)(lattice_step / std::abs(frequency));
    float maximum_error = 0;
    for (int axis = 0; axis < 3; ++axis)
    {
        Coordinate before[] = {x, y, z};
        Coordinate after[] = {x, y, z};
        before[axis] -= h;
        after[axis] += h;
        const float numerical = (noise.GetNoise(after[0], after[1], after[2])
            - noise.GetNoise(before[0], before[1], before[2])) / (2 * h);
        const float error = std::abs(numerical - gradient[axis]);
        maximum_error = std::max(maximum_error, error);
        require(error < 0.0025f * std::abs(frequency) + 0.0002f,
            "analytic gradient disagrees with branch-local numerical derivative");
    }
    return maximum_error;
}

static void source_to_world(const double* q, double* p,
    FastNoiseLite::RotationType3D rotation, float frequency)
{
    // Invert the orthonormal orientation only to place fixtures away from
    // scalar branch switches. Numerical differentiation still uses GetNoise.
    const double a = -0.211324865405187, b = 0.577350269189626;
    const double xy[3][3] = {{1+a,a,b},{a,1+a,b},{-b,-b,b}};
    const double xz[3][3] = {{1+a,b,a},{-b,b,-b},{a,b,1+a}};
    for (int axis = 0; axis < 3; ++axis)
    {
        if (rotation == FastNoiseLite::RotationType3D_None)
            p[axis] = (q[0] + q[1] + q[2]) * (2.0/3.0) - q[axis];
        else
        {
            const double (*matrix)[3] = rotation == FastNoiseLite::RotationType3D_ImproveXYPlanes ? xy : xz;
            p[axis] = matrix[axis][0]*q[0] + matrix[axis][1]*q[1] + matrix[axis][2]*q[2];
        }
        if (frequency != 0) p[axis] /= frequency;
    }
}

int main()
{
    try
    {
        float maximum_error = 0;
        int checked = 0;
        for (int seed : {-1337, 0, 123, 2026})
        {
            FastNoiseLite noise(seed);
            noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2S);
            noise.SetFractalType(FastNoiseLite::FractalType_None);
            for (auto rotation : {FastNoiseLite::RotationType3D_None,
                FastNoiseLite::RotationType3D_ImproveXYPlanes,
                FastNoiseLite::RotationType3D_ImproveXZPlanes})
            {
                noise.SetRotationType3D(rotation);
                for (float frequency : {0.0f, 0.37f, 1.0f, -1.7f})
                {
                    noise.SetFrequency(frequency);
                    for (int i = 0; i < 4; ++i)
                    for (int j = 0; j < 4; ++j)
                    for (int k = 0; k < 4; ++k)
                    {
                        // Odd eighths stay at least 1/8 from lattice/mask
                        // planes. Squared distances to either lattice have an
                        // odd numerator over64, never the support cutoff48/64.
                        // Thus the .001 lattice-space stencil stays inside
                        // one polynomial branch, without selecting fixtures
                        // based on the derivative being tested.
                        const double q[] = {(2*i+1)/8.0-1, (2*j+1)/8.0, (2*k+1)/8.0+1};
                        double p[3];
                        source_to_world(q, p, rotation, frequency);
                        maximum_error = std::max(maximum_error,
                            check_point(noise, p[0], p[1], p[2], frequency));
                        maximum_error = std::max(maximum_error,
                            check_point(noise, (float)p[0], (float)p[1], (float)p[2], frequency));
                        checked += 2;
                    }
                    const auto raw = noise.GetOpenSimplex2SWithGradient(0.13, -0.27, 0.91);
                    noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
                    noise.SetFractalType(FastNoiseLite::FractalType_FBm);
                    const auto independent = noise.GetOpenSimplex2SWithGradient(0.13, -0.27, 0.91);
                    require(raw.value == independent.value && raw.dx == independent.dx
                        && raw.dy == independent.dy && raw.dz == independent.dz,
                        "raw API depends on selected noise/fractal type");
                    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2S);
                    noise.SetFractalType(FastNoiseLite::FractalType_None);
                }
            }
        }
        // Both sides of the inherited a2/a4 selection seam near y=-1.758622.
        // No classical derivative is claimed at the discontinuity itself.
        FastNoiseLite seam(-1337);
        seam.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2S);
        seam.SetFrequency(0.37f);
        for (double y : {-1.76086, -1.75686})
            maximum_error = std::max(maximum_error,
                check_point(seam, -1.02741, y, 0.533843, 0.37f, 0.000037));
        std::cout << "Checked " << checked + 2 << " branch-interior float/double samples; scalar values exact; maximum derivative error "
                  << maximum_error << '\n';
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
