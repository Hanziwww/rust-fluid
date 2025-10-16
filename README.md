# Rust Fluid Simulator

A desktop fluid simulation app built with Rust.

## Run

```bash
# Run in development mode
cargo run

# Run in release mode (best performance)
cargo run --release
```

## Physical Model

This simulator implements a simplified incompressible Navier–Stokes solver:

1. **Advection**: Semi-Lagrangian method
2. **Diffusion**: Ignored (inviscid / Euler flow)
3. **Pressure Projection**: Jacobi iteration to solve the Poisson equation
4. **Boundaries**: Solid walls and fluid boundaries

## Known Limitations

* No multiphase simulation
* Viscosity neglected (best suited for ideal/inviscid flows)

## References

1. **Stam, J.** “Stable Fluids.” *SIGGRAPH 1999 Proceedings.* A seminal unconditionally stable solver popular in graphics; introduces semi-Lagrangian advection with implicit diffusion and projection.
2. **Bridson, R.** *Fluid Simulation for Computer Graphics, 2nd ed.* CRC Press, 2015. Standard textbook covering advection schemes, pressure projection, and practical implementations.
3. **Chorin, A. J.** “Numerical Solution of the Navier–Stokes Equations.” *Math. of Computation,* 1968. Classic fractional-step / projection method used to enforce incompressibility.
4. **Fedkiw, R., Stam, J., Jensen, H. W.** “Visual Simulation of Smoke.” *SIGGRAPH 2001 Proceedings.* Influential graphics paper adding buoyancy and vorticity features atop an Eulerian solver.
5. **Priestley, A.** “A Quasi-Conservative Version of the Semi-Lagrangian Advection Scheme.” *Monthly Weather Review,* 1993. Foundational reference on semi-Lagrangian advection accuracy and conservation in NWP.
