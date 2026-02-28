#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

#include "engine/core/ParticleSoA.hpp"

namespace engine {

struct Vec3 {
    float x{0.0F};
    float y{0.0F};
    float z{0.0F};

    Vec3& operator+=(const Vec3& rhs) noexcept {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }
};

[[nodiscard]] inline Vec3 operator+(Vec3 lhs, const Vec3& rhs) noexcept {
    lhs += rhs;
    return lhs;
}

[[nodiscard]] inline Vec3 operator-(const Vec3& lhs, const Vec3& rhs) noexcept {
    return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

[[nodiscard]] inline Vec3 operator*(const Vec3& lhs, float scalar) noexcept {
    return {lhs.x * scalar, lhs.y * scalar, lhs.z * scalar};
}

class BarnesHutOctree {
public:
    struct Node {
        Vec3 center{};
        float halfSize{0.0F};

        Vec3 centerOfMass{};
        float totalMass{0.0F};

        std::array<int32_t, 8> children{};
        int32_t particleIndex{-1};
        bool isLeaf{true};

        Node() {
            children.fill(-1);
        }
    };

    explicit BarnesHutOctree(float softening = 0.01F) noexcept;

    void build(const ParticleSoA& particles);
    [[nodiscard]] Vec3 accelerationAt(std::size_t index,
                                      const ParticleSoA& particles,
                                      float theta,
                                      float gravitationalConstant) const noexcept;

    [[nodiscard]] const std::vector<Node>& nodes() const noexcept {
        return m_nodes;
    }

private:
    std::vector<Node> m_nodes;
    float m_softening;

    int32_t makeNode(const Vec3& center, float halfSize);
    void insertParticle(int32_t nodeIndex, int32_t particleIndex, const ParticleSoA& particles);
    void subdivide(int32_t nodeIndex);
    [[nodiscard]] int childOctant(const Node& node, Vec3 position) const noexcept;
    [[nodiscard]] Vec3 particlePosition(std::size_t index, const ParticleSoA& particles) const noexcept;
    void computeMassDistribution(int32_t nodeIndex, const ParticleSoA& particles);
};

} // namespace engine
