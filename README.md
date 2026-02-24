# 12 Steps to Navier-Stokes: CFD Solver 

A computational fluid dynamics (CFD) project focused on the numerical solution of the Navier-Stokes equations using Finite Difference Methods (FDM). This repository follows the structured learning path developed by Prof. Lorena A. Barba.

## Project Overview
The goal of this project was to build a CFD solver from scratch, starting from simple 1D linear convection and progressing to the full 2D Navier-Stokes equations for a Cavity Flow and Channel Flow. 

The solver explores key concepts in numerical analysis, such as:
* **Stability Criteria:** Understanding and implementing the Courant-Friedrichs-Lewy (CFL) condition.
* **Discretization:** Applying Forward Time/Backward Space (FTBS) and Central Difference schemes.
* **Pressure-Velocity Coupling:** Solving the Poisson equation for pressure to ensure mass conservation.

## Roadmap & Implementation
The project is divided into 12 incremental steps:

1.  **1D Linear Convection:** The simplest form of the transport equation.
2.  **1D Non-linear Convection:** Introducing the shock wave behavior.
3.  **1D Diffusion:** Implementing the second-order derivative.
4.  **Burgers' Equation:** Combining convection and diffusion.
5.  **2D Linear Convection**
6.  **2D Non-linear Convection**
7.  **2D Diffusion**
8.  **2D Burgers' Equation:** Solving coupled non-linear PDEs.
9.  **2D Laplace Equation:** Foundations for potential flow.
10. **2D Poisson Equation:** Critical for the pressure field.
11. **Cavity Flow with Navier-Stokes:** Modeling a lid-driven cavity.
12. **Channel Flow:** Modeling flow between two parallel plates.

## Technical Stack
* **Language:** Python / MATLAB (Scegli quello che hai usato)
* **Libraries:** NumPy, Matplotlib (Se hai usato Python)
* **Methods:** Finite Difference Method (FDM), Explicit and Implicit schemes.

## Key Results
The final step (Step 11 & 12) simulates the velocity streamlines and pressure contours within a 2D domain, demonstrating the formation of vortices in a lid-driven cavity and the parabolic velocity profile in a channel.

## References
This project is based on the [12 Steps to Navier-Stokes](https://github.com/barbagroup/CFDPython) module by Prof. Lorena A. Barba at Boston University.

---
*Developed as a practical study in Computational Fluid Dynamics.*
