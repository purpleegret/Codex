#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include "engine/core/ParticleSoA.hpp"
#include "engine/physics/gravity/BarnesHutOctree.hpp"

int main() {
    constexpr std::size_t particleCount = 10'000;
    constexpr float G = 6.6743e-11F;

    engine::ParticleSoA particles;
    particles.x.resize(particleCount);
    particles.y.resize(particleCount);
    particles.z.resize(particleCount);
    particles.vx.resize(particleCount, 0.0F);
    particles.vy.resize(particleCount, 0.0F);
    particles.vz.resize(particleCount, 0.0F);
    particles.mass.resize(particleCount, 5.0e20F);

    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-1.0e7F, 1.0e7F);

    for (std::size_t i = 0; i < particleCount; ++i) {
        particles.x[i] = dist(rng);
        particles.y[i] = dist(rng);
        particles.z[i] = dist(rng);
    }

    engine::BarnesHutOctree tree(10.0F);

    const auto buildStart = std::chrono::high_resolution_clock::now();
    tree.build(particles);
    const auto buildEnd = std::chrono::high_resolution_clock::now();

    double checksum = 0.0;
    const auto evalStart = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < particleCount; ++i) {
        const auto a = tree.accelerationAt(i, particles, 0.6F, G);
        checksum += std::sqrt(static_cast<double>(a.x * a.x + a.y * a.y + a.z * a.z));
    }
    const auto evalEnd = std::chrono::high_resolution_clock::now();

    const auto buildMs = std::chrono::duration<double, std::milli>(buildEnd - buildStart).count();
    const auto evalMs = std::chrono::duration<double, std::milli>(evalEnd - evalStart).count();

    std::cout << "Barnes-Hut build time (10k): " << buildMs << " ms\n";
    std::cout << "Acceleration sweep time (10k): " << evalMs << " ms\n";
    std::cout << "Checksum: " << checksum << '\n';

    return 0;
}
