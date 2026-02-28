#pragma once

#include <cstddef>
#include <vector>

namespace engine {

struct ParticleSoA {
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> z;

    std::vector<float> vx;
    std::vector<float> vy;
    std::vector<float> vz;

    std::vector<float> mass;

    [[nodiscard]] std::size_t size() const noexcept {
        return x.size();
    }
};

} // namespace engine
