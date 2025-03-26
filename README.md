# bhevapsolver

This repository contains a Julia implementation of a nondimensionalized version of a first order system of ordinary differential equations describing the evaporation of a charged black hole.

The original set of equations is described in [W.A. Hiscock and L.D. Weems, *Evolution of charged evaporating black holes* Phys.Rev.D 41 (1990) 1142](https://doi.org/10.1103/PhysRevD.41.1142). The nondimensionalized equations implemented here are described in J. Santiago, J. Feng, S. Schuster, and M. Visser, *Immortality through the dark forces: Dark-charge primordial black holes as dark matter candidates*  to appear on the arXiv shortly.

The code requires the 
following Julia packages:

    OrdinaryDiffEq, Plots, Pyplot, ArbNumerics, LaTeXStrings, ImplicitEquations, Serialization

To run the code, open a terminal and cd to the following directory:

    > cd src/scripts/

To run the HWscript_DE.jl script, run the following command:

    > julia HWscript_DE.jl

Alternatively, you can simply paste the contents of HWscript_DE.jl into the julia REPL (make sure julia is running in the "src/scripts/" directory).

