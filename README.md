# MixingIBAMR

Diffusive mixing in binary fluid mixtures in the presence of thermal fluctuations
Aleksandar Donev (donev@courant.nyu.edu), Courant Institute of Mathematical Sciences
and Boyce Griffith (boyceg@gmail.com), University of North Carolina, Chapel Hill

This software relies on the IBAMR library (https://github.com/IBAMR) for the fluid dynamics solvers.

These codes implement numerical methods for simulating giant nonequilibrium fluctuations 
in binary fluid mixtures, as described in a sequence of papers. The primary source is 
the description of the temporal integrators and tests described in (see doc/MultiscaleIntegrators.pdf)

[1] "Multiscale temporal integrators for fluctuating hydrodynamics"
S. Delong, Y. Sun, B.E. Griffith, E. Vanden-Eijnden and A. Donev
Phys. Rev. E, 90, 063312, 2014, http://arxiv.org/abs/1410.0240.
(note: please use the doc or ArXiv file instead of PRE since some equation typos have been fixed)

In particular, these codes solve either the:
1) Inertial equations (63,64) in this paper, or
2) Overdamped (inertia-free) equations (65)

Note: Make sure you are using the correct main*.C (inertial or overdamped) for each input file (see below).
The temporal integrators here are modified from the basic IBAMR library because we often need to swap
the order in which velocity and concentration are updated. This is done by partially reusing the exisiting
temporal integrators in IBAMR, and partially updating them. The result is that setting the correct options
and understanding their interactions is quite complex even for ourselves.
If you want to use or understand this code make sure to read doc/README carefully.

The spatial discretization, and in particular details about the handling of the stochastic forcing
in the fluid equation, can be found in:

[2] "Staggered Schemes for Fluctuating Hydrodynamics"
F. Balboa and J. Bell and R. Delgado-Buscalioni and A. Donev and T. G. Fai and B. Griffith and C. Peskin
SIAM J. Multiscale Modeling and Simulation, 10(4):1369-1408, 2012, http://arxiv.org/abs/1108.5188.
(but note that this code only implements the incompressible equations, see fluam code for compressible)

--------------------
A. The code main.C

Is meant to provide integrators for both overdamped and inertial in the case when there is no gravity, i.e.,
when the concentration/temperature (scalar) equation does not feed back into the velocity equation
In this case we can avoid a second fluid solve and gain efficiency while still being second-order accurate.
This code also provides access to the default IBAMR integrators for advection-diffusion, which
can be second-order accurate also in the presence of both inertia and gravity.
For overdamped limit with gravity see main_PC.C below.

i) In the overdamped case, in the absence of gravity,
set use_split_time_stepping=TRUE and CENTER_U_ADV_DIFF=FALSE.
This implements a version of Algorithm 3 in paper [1] streamlined for the case g=0
by avoiding the second Stokes solve, which is identical to the first.
Note that the Soret flux is handled implicitly in the concentration solver by modifying the Helmholtz solver.
More physical details as well as comparisons with actual experiments can be found in:

[3] "Dynamic scaling for the growth of non-equilibrium fluctuations during thermophoretic diffusion in microgravity"
R. Cerbino, Y. Sun, A. Donev and A. Vailati
Submitted, 2015, http://arxiv.org/abs/1502.03693.

See examples/GRADFLEX-transient-?d.input for dynamical evolution of fluctuations.
For steady-state spectrum calculations (static structure factor), see
examples/GRADFLEX-dynamic-3d.input

ii) For inertial, if there is no gravity, set CENTER_U_ADV_DIFF=TRUE, which
gives an algorithm similar to that proposed in paper [2] but simplified for g=0.
See examples/GRADFLEX-temperature-2d.input for temperature equation in paper [3].

iii) For inertial, if there is gravity (the same of course also works if g=0), 
set use_split_time_stepping=FALSE to obtain Algorithm 1 in paper [1],
as detailed in Section V.B in paper [1] above, see examples/THN-C12-20K-inertial.input.
More physical details as well as comparisons with actual experiments can be found in:

[4] "Slowing-down of non-equilibrium concentration fluctuations in confinement"
C. Giraudet and H. Bataller and Y. Sun and A. Donev and J. M. O. de Zárate and F. Croccolo
Submitted, 2015, http://arxiv.org/abs/1410.6524.

--------------------
B. Overdamped code mainPC.C

Is meant to provide second-order accuracy also for the overdamped limit in the presence of gravity.
It solves the overdamped equations using Algorithm 3 in paper [1].
See examples/THN-C12-20K-overdamped.input for a comparison with inertial, as discussed in paper [4].
