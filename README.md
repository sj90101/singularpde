# singularpde

This repo contains Matlab functions for solving some partial differential equations and integral equations whose solutions are singular or nearly singular.
The singularity of the solution is induced by the geometry and/or the right-hand-side functions. The method used is the RCIP (Recursively compressed inverse preconditioning) 
method. This is joint work with Johan Helsing at Lund University in Sweden.

# Callable Matlab functions

* inireg3j.m, modex3j.m, modgp1j.m: solves the first Dirichlet problem of the biharmonic and modified biharmonic equations in domains with corners in two dimensions.
* bgkw12j.m: solves the linearized BGKW equation for steady Couette flow.
* testhelmos.m: solves the Helmholtz Dirichlet and Neumann problems on piecewise smooth open curves.

# Citing

If you find these codes useful in your work, please star this repository and cite it and one of the following papers. 


```
@article {HJ2018modbihar,
    AUTHOR = {Helsing, Johan and Jiang, Shidong},
     TITLE = {On integral equation methods for the first {D}irichlet problem
              of the biharmonic and modified biharmonic equations in
              nonsmooth domains},
  JOURNAL = {SIAM Journal on Scientific Computing},
    VOLUME = {40},
      YEAR = {2018},
    NUMBER = {4},
     PAGES = {A2609--A2630},
       DOI = {10.1137/17M1162238},
       URL = {https://doi.org/10.1137/17M1162238},
}
```
```
@article{HJ2022singularrhs,
title = {Solving {F}redholm second-kind integral equations with singular right-hand sides on non-smooth boundaries},
author = {Johan Helsing and Shidong Jiang},
journal = {Journal of Computational Physics},
volume = {448},
pages = {110714},
year = {2022},
issn = {0021-9991},
doi = {https://doi.org/10.1016/j.jcp.2021.110714},
url = {https://www.sciencedirect.com/science/article/pii/S0021999121006094},
}
```
```
@misc{HJ2024helmos,
      title={The {H}elmholtz {D}irichlet and {N}eumann problems on piecewise smooth open curves}, 
      author={Johan Helsing and Shidong Jiang},
      year={2024},
      eprint={2411.05761},
      archivePrefix={arXiv},
      primaryClass={math.NA},
      note={\url{https://arxiv.org/abs/2411.05761}},
}
```

# Main developers

* Johan Helsing, Lund University, Sweden
* Shidong Jiang, Flatiron Institute, Simons Foundation

