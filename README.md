# MixingIBAMR

Diffusive mixing in binary fluid mixtures in the presence of thermal fluctuations
Aleksandar Donev (donev@courant.nyu.edu), Courant Institute of Mathematical Sciences
and Boyce Griffith (boyceg@gmail.com), University of North Carolina, Chapel Hill

This software relies on the IBAMR library (https://github.com/IBAMR) for the fluid dynamics solvers.

These codes implement numerical methods for simulating giant nonequilibrium fluctuations 
in binary fluid mixtures, as described in a sequence of papers. The primary source is 
the description of the temporal integrators and tests described in (see doc/MultiscaleIntegrators.pdf)

1. "Multiscale temporal integrators for fluctuating hydrodynamics"
S. Delong, Y. Sun, B.E. Griffith, E. Vanden-Eijnden and A. Donev
Phys. Rev. E, 90, 063312, 2014, http://arxiv.org/abs/1410.0240.
(note: please use the doc or ArXiv file instead of PRE since some equation typos have been fixed)

The spatial discretization, and in particular details about the handling of the stochastic forcing
in the fluid equation, can be found in:

2. "Staggered Schemes for Fluctuating Hydrodynamics"
F. Balboa and J. Bell and R. Delgado-Buscalioni and A. Donev and T. G. Fai and B. Griffith and C. Peskin
SIAM J. Multiscale Modeling and Simulation, 10(4):1369-1408, 2012 [arXiv:1108.5188]
(but note that this code only implements the incompressible equations, see fluam code for compressible)

--------------------
A. Inertial code main.C

solves the inertial equations (63,64) using Algorithm 1,
as detailed in Section V.B in paper [1] above, see examples/EARTH-SORET-dynamic-2d.input
More physical details as well as comparisons with actual experiments can be found in:

3. "Slowing-down of non-equilibrium concentration fluctuations in confinement"
C. Giraudet and H. Bataller and Y. Sun and A. Donev and J. M. O. de Zárate and F. Croccolo
Submitted, 2015 [ArXiv:1410.6524]

see examples/THN-C12-20K-inertial.input for inertial
and examples/THN-C12-20K-overdamped.input for a comparison with overdamped

Note that part of the temporal integrator control is provided via the IBAMR input file,
but part of it is overwritten manually in main.C.

--------------------
B. Overdamped code mainPC.C

solves the overdamped (inertia-free) equations (65) using Algorithm 3,
as detailed in Section V.A in paper [1] above, see examples/GRADFLEX-transient-?d.input
More physical details as well as comparisons with actual experiments can be found in:

"Dynamic scaling for the growth of non-equilibrium fluctuations during thermophoretic diffusion in microgravity"
R. Cerbino, Y. Sun, A. Donev and A. Vailati
Submitted, 2015 [ArXiv:1502.03693].

See examples/GRADFLEX-experiment-3d.input
For steady-state spectrum calculations (static structure factor), see
examples/GRADFLEX-dynamic-3d.input for concentration
examples/GRADFLEX-temperature-2d.input for temperature

Note that part of the temporal integrator control is provided via the IBAMR input file,
but part of it is overwritten manually in mainPC.C.
Also note that the Soret flux is handled implicitly in the concentration solver by modifying the Helmholtz solver.
