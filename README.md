# Codex - Celestial SPH/N-Body Prototype

This repository now contains a high-performance-oriented 3D Barnes-Hut octree implementation and data-oriented particle storage as the first milestone of a full SPH + N-body celestial simulator.

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

## Run benchmark demo

```bash
./build/octree_demo
```

The demo creates 10,000 particles, builds the octree, evaluates accelerations for all particles, and prints timing/checksum values.

See `docs/architecture.md` for the full target project structure and physics foundations (SPH kernels, EOS, artificial viscosity, leapfrog integration).
