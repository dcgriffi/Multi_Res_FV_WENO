# Multi_Res_FV_WENO

This repository contains a MATLAB implementation of a finite volume (FV) multi-resolution Weighted Essentially Non-Oscillatory (WENO) scheme for hyperbolic conservation laws, applied on cell-vertex grids. The code supports 1D and 2D domains, including rectangular and cell-vertex mesh structures.

## Overview

The project focuses on solving hyperbolic conservation laws using a multi-resolution approach combined with constrained least squares (LS) reconstruction techniques. It demonstrates how oscillations form near discontinuities or shocks and how they can be controlled effectively.

## File Structure and Description

- **LS_*.m**  
  These scripts implement constrained least squares reconstructions—both quadratic and linear—for scalar functions. They serve as the building blocks for the multi-resolution algorithm.

- **FV_LS_*.m**  
  These contain the quadratic least squares reconstructions applied within a finite volume framework. They illustrate the formation of oscillations near discontinuous initial conditions or shock regions.

- **Multires_*.m / .asv**  
  These files implement the multi-resolution algorithm for hyperbolic conservation laws:  
  - *Multires_1D.m*: 1D cases for general flux function \( f \).  
  - *Multires_2D_Rectangle.m*: 2D rectangular grid cases.  
  - *Multires_2D_Cell_Vertex.m*: 2D cell-vertex grid cases.  
  - *Multires_2D_Cell_Vertex_adv.m*: Advection-specific 2D cell-vertex test cases showing convergence and oscillation control. This is the primary finalized output of the project.

- **Mesh and Utility Functions:**  
  - `gen_reg_tri_mesh.m`: Generates regular triangular meshes.  
  - `get_fv_arrays.m`: Prepares finite volume grid arrays and connectivity.  
  - `integQuad.m`: Performs quadrature integrations.  
  - `lgwt.m`: Computes Legendre-Gauss nodes and weights.

- **Additional Examples:**  
  - `solveBurgers2.m`: Solves Burgers' equation as a test case.  
  - `LS_Lin_2D_Rectangle.m` & `LS_Quad_2D_Rectangle.m`: Least squares for 2D rectangular grids.

## Features

- Multi-resolution WENO scheme for hyperbolic conservation laws.  
- Supports cell-vertex and rectangular grids in 1D and 2D.  
- Constrained least squares reconstruction techniques (linear and quadratic).  
- Demonstrates oscillation formation and control near shocks.  
- Provides advection-specific convergence tests.

## Usage

To run a test case, start with one of the *Multires* scripts (e.g., `Multires_2D_Cell_Vertex_adv.m`) to observe the solver behavior on advection problems. Modify the parameters and initial conditions inside these scripts to explore different scenarios.

## Requirements

- MATLAB R2018b or later (for better compatibility with function syntax).  
- No external toolboxes required.
