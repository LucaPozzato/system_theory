# System Theory Project

This repository contains MATLAB code and LaTeX documentation for the System Theory project. The project involves the analysis of a nonlinear dynamic system, including stability, bifurcations, chaos, and Lyapunov exponents.

## Project Structure

### MATLAB Code

The MATLAB code is organized into several scripts and functions, each focusing on different aspects of the system analysis:

- `matlab/sys_bifurcation.m`: Analyzes the bifurcations of the system by computing equilibria and their stability over a range of parameters.
- `matlab/sys.m`: Solves for the equilibria of the system and evaluates their stability using the Jacobian matrix.
- `matlab/sys_lyapunov.m`: Computes the Lyapunov exponents of the system to determine its chaotic behavior.
- `matlab/sys_stability.m`: Evaluates the stability of the system's equilibria and simulates limit cycles.
- `matlab/sys_chaos.m`: Simulates the system with time-varying resistance to analyze chaotic behavior and generates a Poincaré map.
- `matlab/sys_div.m`: Computes and visualizes the divergence of the system for different parameter values.
- `matlab/sys_animation.m`: Creates an animation of the system's behavior over time with two slightly different initial conditions.

### Documentation

The pdf document `documents/document.pdf` provides a detailed analysis of the system, including:

- State equations derived from the electrical circuit.
- Analysis of oscillations with large resistance values.
- Stability and bifurcation analysis of the equilibria.
- Invariants of the system as parameters vary.
- Number of bifurcations and their types.
- Analysis of chaos, Poincaré sections, and Lyapunov exponents.
- Bifurcation diagrams and qualitative sketches of the system's behavior.
