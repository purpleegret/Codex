# 3D Celestial SPH/N-Body Simulator Architecture

## 1) Modular project file structure

```text
CelestialSPHSim/
├── CMakeLists.txt
├── README.md
├── docs/
│   └── architecture.md
├── include/
│   └── engine/
│       ├── core/
│       │   └── ParticleSoA.hpp              # Data-oriented particle storage (SoA)
│       ├── physics/
│       │   ├── gravity/
│       │   │   └── BarnesHutOctree.hpp      # O(n log n) gravity accelerator
│       │   └── sph/
│       │       ├── CubicSplineKernel.hpp    # W(r,h), ∇W for SPH density/pressure
│       │       ├── EquationOfState.hpp      # Material EOS (e.g., Tait/Murnaghan)
│       │       └── ArtificialViscosity.hpp  # Collision shock handling
│       ├── integrators/
│       │   └── LeapfrogIntegrator.hpp       # Symplectic kick-drift-kick stepping
│       └── render/
│           ├── ParticleRenderer.hpp         # Instanced rendering path
│           └── GpuBuffers.hpp               # GPU upload and compute interop
└── src/
    ├── app/
    │   └── main.cpp                         # Benchmark/demo entry point
    ├── physics/
    │   ├── gravity/
    │   │   └── BarnesHutOctree.cpp          # Octree build + force traversal
    │   └── sph/
    │       ├── CubicSplineKernel.cpp
    │       ├── EquationOfState.cpp
    │       └── ArtificialViscosity.cpp
    ├── integrators/
    │   └── LeapfrogIntegrator.cpp
    └── render/
        └── ParticleRenderer.cpp
```

> This repository currently implements the gravity octree core and SoA particle layout, while the SPH/integrator/renderer modules are scaffold targets for the next tasks.

## 2) Mathematical foundation

### 2.1 Barnes-Hut gravity (3D octree)

- Build a cubic root cell that bounds all particles.
- Recursively subdivide into 8 octants until each leaf holds at most one particle.
- For each node, compute aggregate mass and center of mass:

\[
M_n = \sum_{i \in n} m_i, \qquad
\mathbf{r}_{\text{cm},n} = \frac{1}{M_n} \sum_{i \in n} m_i\,\mathbf{r}_i
\]

- Use the opening criterion during traversal:

\[
\frac{s}{r} < \theta
\]

where \(s\) is node size, \(r\) is distance from query particle to node COM, and \(\theta\) is accuracy parameter.

- If criterion passes, approximate the node as one pseudo-body:

\[
\mathbf{a}_i = G M_n \frac{\mathbf{r}_{\text{cm},n} - \mathbf{r}_i}{\left(\|\mathbf{r}_{\text{cm},n} - \mathbf{r}_i\|^2 + \varepsilon^2\right)^{3/2}}
\]

with softening \(\varepsilon\) to avoid singular acceleration.

### 2.2 SPH with cubic spline kernel

Density at particle \(i\):
\[
\rho_i = \sum_j m_j W(\|\mathbf{r}_i-\mathbf{r}_j\|, h)
\]

Pressure force (symmetric form):
\[
\frac{d\mathbf{v}_i}{dt}\Big|_{p} = -\sum_j m_j\left(\frac{P_i}{\rho_i^2}+\frac{P_j}{\rho_j^2}\right)\nabla W_{ij}
\]

Cubic spline kernel in 3D with \(q=r/h\):
\[
W(r,h)=\frac{1}{\pi h^3}
\begin{cases}
1-\frac{3}{2}q^2+\frac{3}{4}q^3,&0\le q<1\\
\frac{1}{4}(2-q)^3,&1\le q<2\\
0,&q\ge2
\end{cases}
\]

### 2.3 Equation of state (material compressibility/strength)

A common planetary-fluid starting point is Tait-like EOS:
\[
P = B\left[\left(\frac{\rho}{\rho_0}\right)^\gamma - 1\right]
\]

- \(\rho_0\): reference density
- \(B\): bulk modulus-like stiffness
- \(\gamma\): compressibility exponent

For rocky/metallic regimes, this can later be replaced with Murnaghan/Tillotson/ANEOS tables for impact physics.

### 2.4 Artificial viscosity (shock/collision handling)

Monaghan-type artificial viscosity stabilizes shocks:
\[
\Pi_{ij}=
\begin{cases}
\frac{-\alpha c_{ij}\mu_{ij}+\beta \mu_{ij}^2}{\bar{\rho}_{ij}}, & \mathbf{v}_{ij}\cdot\mathbf{r}_{ij}<0\\
0,&\text{otherwise}
\end{cases}
\]

\[
\mu_{ij} = \frac{h\,\mathbf{v}_{ij}\cdot\mathbf{r}_{ij}}{\|\mathbf{r}_{ij}\|^2+\eta^2}
\]

It is added to momentum and energy equations to capture impact dissipation without explicit Riemann solvers.

### 2.5 Symplectic leapfrog integration

Kick-drift-kick:
1. \(\mathbf{v}^{n+1/2}=\mathbf{v}^n + \frac{\Delta t}{2}\mathbf{a}(\mathbf{x}^n)\)
2. \(\mathbf{x}^{n+1}=\mathbf{x}^n + \Delta t\,\mathbf{v}^{n+1/2}\)
3. Recompute \(\mathbf{a}(\mathbf{x}^{n+1})\)
4. \(\mathbf{v}^{n+1}=\mathbf{v}^{n+1/2} + \frac{\Delta t}{2}\mathbf{a}(\mathbf{x}^{n+1})\)

This symplectic form preserves long-term orbital behavior better than explicit Euler or RK for Hamiltonian systems.

## 3) Data-oriented design notes

- Structure-of-arrays (SoA) storage for particle attributes maximizes cache line utilization and SIMD potential.
- Octree nodes are stored in a contiguous node pool (`std::vector<Node>`) to reduce pointer chasing.
- Traversal is iterative for force evaluation (explicit stack), avoiding recursion overhead in hot loops.
- The same SoA layout maps naturally to GPU SSBO/UAV buffers for compute-shader acceleration.
