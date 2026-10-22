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
            "analytic gradient disagrees with numerical derivative");
    }
    return maximum_error;
}

static void source_to_world(const double* q, double* p,
    FastNoiseLite::RotationType3D rotation, float frequency)
{
    // Invert the sampler's orthonormal orientation to construct source-space
    // fixtures; derivatives are still checked through the public API.
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

static const double GradientScale = 0.75 * 0.75 * 0.75 * 0.75 * 9.046026385208288;

static double dot3(const double* a, const double* b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void source_gradient(const FastNoiseLite::NoiseGradient3D& sample, double* g)
{
    // The oracle uses frequency one and the symmetric default orientation.
    const double r = (sample.dx + sample.dy + sample.dz) * (2.0 / 3.0);
    g[0] = (r - sample.dx) / GradientScale;
    g[1] = (r - sample.dy) / GradientScale;
    g[2] = (r - sample.dz) / GradientScale;
}

static float complete_support_value(const double* q, int seed)
{
    FastNoiseLite nodes(seed);
    nodes.SetFrequency(1.0f);
    double sum=0;
    for (int lattice=0;lattice<2;++lattice) {
        const double offset=lattice ? 0.5 : 0.0;
        const int base[3] = {(int)std::floor(q[0]+offset),(int)std::floor(q[1]+offset),(int)std::floor(q[2]+offset)};
        for (int x=0;x<2;++x) for (int y=0;y<2;++y) for (int z=0;z<2;++z) {
            const double node[3] = {base[0]+x-offset,base[1]+y-offset,base[2]+z-offset};
            const double d[3] = {q[0]-node[0],q[1]-node[1],q[2]-node[2]};
            const double attenuation=0.75-dot3(d,d);
            if (attenuation<=0) continue;
            double world_node[3];
            source_to_world(node, world_node, FastNoiseLite::RotationType3D_None, 1);
            const auto sample = nodes.GetOpenSimplex2SWithGradient(
                world_node[0], world_node[1], world_node[2]);
            double g[3];
            source_gradient(sample, g);
            require(std::abs(dot3(g,g)-2.0)<0.0002,"recovered nodal gradient has unexpected norm");
            sum += 9.046026385208288*attenuation*attenuation*attenuation*attenuation*dot3(g,d);
        }
    }
    return (float)sum;
}

int main()
{
    try
    {
        float maximum_error = 0;
        int checked = 0;
        int oracle_checked = 0;
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
        // Transport the known old source seam into every orientation, rather
        // than treating its old world coordinates as source coordinates.
        const double old_point[3] = {-1.02741, -1.75862226094, 0.533843};
        double seam_source[3];
        source_to_world(old_point, seam_source, FastNoiseLite::RotationType3D_None, 1);
        for (int axis = 0; axis < 3; ++axis) seam_source[axis] *= 0.37f;
        const double seam_step[3] = {2.0/3 * 0.37f, -1.0/3 * 0.37f, 2.0/3 * 0.37f};
        int seam_checked = 0;
        for (int seed : {-1337, 0, 123, 2026})
        {
            FastNoiseLite seam(seed);
            seam.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2S);
            seam.SetFrequency(0.37f);
            for (auto rotation : {FastNoiseLite::RotationType3D_None,
                FastNoiseLite::RotationType3D_ImproveXYPlanes,
                FastNoiseLite::RotationType3D_ImproveXZPlanes})
            {
                seam.SetRotationType3D(rotation);
                for (int permutation = 0; permutation < 3; ++permutation)
                {
                    double q[3], lower[3], upper[3], p[3], before[3], after[3];
                    for (int axis = 0; axis < 3; ++axis)
                    {
                        const int index = (axis + permutation) % 3;
                        q[axis] = seam_source[index];
                        lower[axis] = q[axis] - 1e-6 * seam_step[index];
                        upper[axis] = q[axis] + 1e-6 * seam_step[index];
                    }
                    source_to_world(q, p, rotation, 0.37f);
                    source_to_world(lower, before, rotation, 0.37f);
                    source_to_world(upper, after, rotation, 0.37f);
                    const auto center = seam.GetOpenSimplex2SWithGradient(p[0], p[1], p[2]);
                    const auto lo = seam.GetOpenSimplex2SWithGradient(before[0], before[1], before[2]);
                    const auto hi = seam.GetOpenSimplex2SWithGradient(after[0], after[1], after[2]);
                    require(std::abs(hi.value - lo.value) < 2e-5, "former seam scalar is discontinuous");
                    for (auto neighbor : {lo, hi})
                        require(std::abs(neighbor.dx - center.dx) < 0.0015
                            && std::abs(neighbor.dy - center.dy) < 0.0015
                            && std::abs(neighbor.dz - center.dz) < 0.0015,
                            "former seam gradient is discontinuous");
                    for (const double* point : {before, p, after})
                        maximum_error = std::max(maximum_error,
                            check_point(seam, point[0], point[1], point[2], 0.37f));
                    ++seam_checked;
                }
            }
        }
        FastNoiseLite oracle(1337);
        oracle.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2S);
        oracle.SetRotationType3D(FastNoiseLite::RotationType3D_None);
        oracle.SetFrequency(1.0f);
        for (int i=0;i<160;++i) {
            double q[3] = {std::sin(i*1.31)*4.7, std::cos(i*0.73)*3.9-2.0, std::sin(i*0.37)*5.1};
            if (i%4==0) q[i%3]=std::floor(q[i%3])+0.5;
            if (i%4==1) q[i%3]=std::floor(q[i%3]);
            double p[3]; source_to_world(q,p,FastNoiseLite::RotationType3D_None,1);
            const float actual=oracle.GetNoise(p[0],p[1],p[2]);
            const float expected=complete_support_value(q,1337);
            require(std::abs(actual-expected)<5e-6f,"public scalar disagrees with complete 16-vertex support oracle");
            ++oracle_checked;
        }
        for (double transition : {0.0, 0.5, std::sqrt(0.75)-1e-7,
            std::sqrt(0.75), std::sqrt(0.75)+1e-7}) {
            const double q[3]={transition,0.0,0.0};
            double p[3]; source_to_world(q,p,FastNoiseLite::RotationType3D_None,1);
            const float actual=oracle.GetNoise(p[0],p[1],p[2]);
            const float expected=complete_support_value(q,1337);
            require(std::abs(actual-expected)<5e-6f,"support/tie fixture disagrees with complete support oracle");
            ++oracle_checked;
        }
        std::cout << "Checked " << checked << " float/double samples, " << seam_checked
                  << " former-seam fixtures, and " << oracle_checked
                  << " independent complete-support points; fused/value-only values exact; maximum derivative error "
                  << maximum_error << '\n';
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
