#include "engine/physics/gravity/BarnesHutOctree.hpp"

#include <algorithm>
#include <cmath>
#include <stack>

namespace engine {

BarnesHutOctree::BarnesHutOctree(float softening) noexcept
    : m_softening(softening) {
}

int32_t BarnesHutOctree::makeNode(const Vec3& center, float halfSize) {
    Node node;
    node.center = center;
    node.halfSize = halfSize;
    m_nodes.push_back(node);
    return static_cast<int32_t>(m_nodes.size() - 1);
}

Vec3 BarnesHutOctree::particlePosition(std::size_t index, const ParticleSoA& particles) const noexcept {
    return {particles.x[index], particles.y[index], particles.z[index]};
}

int BarnesHutOctree::childOctant(const Node& node, Vec3 position) const noexcept {
    int octant = 0;
    octant |= static_cast<int>(position.x >= node.center.x) << 0;
    octant |= static_cast<int>(position.y >= node.center.y) << 1;
    octant |= static_cast<int>(position.z >= node.center.z) << 2;
    return octant;
}

void BarnesHutOctree::subdivide(int32_t nodeIndex) {
    const Vec3 parentCenter = m_nodes[nodeIndex].center;
    const float childHalf = m_nodes[nodeIndex].halfSize * 0.5F;

    for (int i = 0; i < 8; ++i) {
        const float sx = (i & 1) ? 1.0F : -1.0F;
        const float sy = (i & 2) ? 1.0F : -1.0F;
        const float sz = (i & 4) ? 1.0F : -1.0F;
        const Vec3 childCenter{
            parentCenter.x + sx * childHalf,
            parentCenter.y + sy * childHalf,
            parentCenter.z + sz * childHalf,
        };
        m_nodes[nodeIndex].children[i] = makeNode(childCenter, childHalf);
    }

    m_nodes[nodeIndex].isLeaf = false;
}

void BarnesHutOctree::insertParticle(int32_t nodeIndex,
                                     int32_t particleIndex,
                                     const ParticleSoA& particles) {
    if (m_nodes[nodeIndex].isLeaf && m_nodes[nodeIndex].particleIndex == -1) {
        m_nodes[nodeIndex].particleIndex = particleIndex;
        return;
    }

    if (m_nodes[nodeIndex].isLeaf) {
        const int32_t existing = m_nodes[nodeIndex].particleIndex;
        m_nodes[nodeIndex].particleIndex = -1;
        subdivide(nodeIndex);

        const Node parent = m_nodes[nodeIndex];
        const Vec3 existingPos = particlePosition(existing, particles);
        const int existingOct = childOctant(parent, existingPos);
        insertParticle(parent.children[existingOct], existing, particles);
    }

    const Node parent = m_nodes[nodeIndex];
    const Vec3 pos = particlePosition(particleIndex, particles);
    const int oct = childOctant(parent, pos);
    insertParticle(parent.children[oct], particleIndex, particles);
}

void BarnesHutOctree::computeMassDistribution(int32_t nodeIndex, const ParticleSoA& particles) {
    Node& node = m_nodes[nodeIndex];
    node.totalMass = 0.0F;
    node.centerOfMass = {};

    if (node.isLeaf) {
        if (node.particleIndex >= 0) {
            const auto i = static_cast<std::size_t>(node.particleIndex);
            node.totalMass = particles.mass[i];
            node.centerOfMass = particlePosition(i, particles);
        }
        return;
    }

    for (const int32_t childIndex : node.children) {
        if (childIndex < 0) {
            continue;
        }

        computeMassDistribution(childIndex, particles);
        const Node& child = m_nodes[childIndex];
        node.totalMass += child.totalMass;
        node.centerOfMass += child.centerOfMass * child.totalMass;
    }

    if (node.totalMass > 0.0F) {
        const float invMass = 1.0F / node.totalMass;
        node.centerOfMass = node.centerOfMass * invMass;
    }
}

void BarnesHutOctree::build(const ParticleSoA& particles) {
    m_nodes.clear();
    if (particles.size() == 0) {
        return;
    }

    m_nodes.reserve(particles.size() * 2);

    const auto [minXIt, maxXIt] = std::minmax_element(particles.x.begin(), particles.x.end());
    const auto [minYIt, maxYIt] = std::minmax_element(particles.y.begin(), particles.y.end());
    const auto [minZIt, maxZIt] = std::minmax_element(particles.z.begin(), particles.z.end());

    const Vec3 minCorner{*minXIt, *minYIt, *minZIt};
    const Vec3 maxCorner{*maxXIt, *maxYIt, *maxZIt};

    const Vec3 center{
        0.5F * (minCorner.x + maxCorner.x),
        0.5F * (minCorner.y + maxCorner.y),
        0.5F * (minCorner.z + maxCorner.z),
    };

    const float spanX = maxCorner.x - minCorner.x;
    const float spanY = maxCorner.y - minCorner.y;
    const float spanZ = maxCorner.z - minCorner.z;
    const float halfSize = 0.5F * std::max({spanX, spanY, spanZ}) + std::numeric_limits<float>::epsilon();

    const int32_t root = makeNode(center, halfSize);

    for (std::size_t i = 0; i < particles.size(); ++i) {
        insertParticle(root, static_cast<int32_t>(i), particles);
    }

    computeMassDistribution(root, particles);
}

Vec3 BarnesHutOctree::accelerationAt(std::size_t index,
                                     const ParticleSoA& particles,
                                     float theta,
                                     float gravitationalConstant) const noexcept {
    if (m_nodes.empty()) {
        return {};
    }

    const Vec3 query = particlePosition(index, particles);
    Vec3 acc{};

    std::stack<int32_t> stack;
    stack.push(0);

    while (!stack.empty()) {
        const int32_t nodeIndex = stack.top();
        stack.pop();

        const Node& node = m_nodes[nodeIndex];
        if (node.totalMass <= 0.0F) {
            continue;
        }

        if (node.isLeaf && node.particleIndex == static_cast<int32_t>(index)) {
            continue;
        }

        const Vec3 dr = node.centerOfMass - query;
        const float r2 = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z + m_softening * m_softening;
        const float r = std::sqrt(r2);
        const float sizeOverDist = (2.0F * node.halfSize) / r;

        if (node.isLeaf || sizeOverDist < theta) {
            const float invR3 = 1.0F / (r2 * r);
            const float scale = gravitationalConstant * node.totalMass * invR3;
            acc += dr * scale;
            continue;
        }

        for (const int32_t childIndex : node.children) {
            if (childIndex >= 0) {
                stack.push(childIndex);
            }
        }
    }

    return acc;
}

} // namespace engine
