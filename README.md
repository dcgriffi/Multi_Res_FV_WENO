# Multi_Res_FV_WENO

This repository contains a MATLAB implementation of a finite volume (FV) multi-resolution Weighted Essentially Non-Oscillatory (WENO) scheme for hyperbolic conservation laws, applied on cell-vertex grids. The code supports 1D and 2D domains, including rectangular and cell-vertex mesh structures.

## Background

This code was developed to explore multi-resolution WENO reconstruction as a less expensive alternative to flux-corrected transport (FCT) for tracer advection in the FESOM2 ocean model [1], which uses a median-dual cell-vertex discretization. The multi-resolution WENO schemes of Zhu and Shu [2, 3] use nested centered stencils to obtain high-order reconstruction in smooth regions and first-order reconstruction near discontinuities, which controls spurious oscillations. This code extends that approach to cell-vertex meshes, removes the intermediate linear stage of the original algorithm, and tests alternative smoothness indicators.

The code was written by D. Chloe Griffin. The project and its underlying theory were developed in collaboration with Chi-Wang Shu (Brown University) and Sergey Danilov (Alfred Wegener Institute) for the Alfred Wegener Institute.

## Overview

The project focuses on solving hyperbolic conservation laws using a multi-resolution approach combined with constrained least squares (LS) reconstruction techniques. It demonstrates how oscillations form near discontinuities or shocks and how they can be controlled effectively.

## File Structure and Description

- **LS_*.m**  
  These scripts implement constrained least squares reconstructions—both quadratic and linear—for scalar functions. They serve as the building blocks for the multi-resolution algorithm.

- **FV_LS_*.m**  
  These contain the quadratic least squares reconstructions applied within a finite volume framework. They illustrate the formation of oscillations near discontinuous initial conditions or shock regions and third order convergence in smooth regions.

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
  - `solveBurgers2.m`: Solves Burgers' equation as a test case.  

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

## References

1. S. Danilov, D. Sidorenko, Q. Wang, T. Jung, *The finite-volume sea ice–ocean model (FESOM2)*, Geoscientific Model Development **10** (2017) 765–789. https://doi.org/10.5194/gmd-10-765-2017
2. J. Zhu, C.-W. Shu, *A new type of multi-resolution WENO schemes with increasingly higher order of accuracy*, Journal of Computational Physics **375** (2018) 659–683. https://doi.org/10.1016/j.jcp.2018.09.003
3. J. Zhu, C.-W. Shu, *A new type of multi-resolution WENO schemes with increasingly higher order of accuracy on triangular meshes*, Journal of Computational Physics **392** (2019) 19–33. https://doi.org/10.1016/j.jcp.2019.04.027

## License

This code is released under the MIT License; see [LICENSE](LICENSE).

Third-party code included in this repository keeps its original authorship:

- `lgwt.m` was written by Greg von Winckel (2004) and is distributed via MATLAB
  File Exchange.
- The mesh and stencil routines (`gen_reg_tri_mesh.m`, `get_fv_arrays.m`, and
  `get_mesh_sten_info` in `Multires_2D_Cell_Vertex_adv.m`) are modified from code
  by Sergey Danilov and are shared with his permission.
