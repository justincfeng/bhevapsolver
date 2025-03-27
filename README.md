# bhevapsolver

This repository contains a Julia language implementation of a nondimensionalized version of a first order system of ordinary differential equations describing the evaporation of a charged black hole.

The original set of equations may be found in [W.A. Hiscock and L.D. Weems, *Evolution of charged evaporating black holes* Phys.Rev.D 41 (1990) 1142](https://doi.org/10.1103/PhysRevD.41.1142), and the nondimensionalized equations implemented here are described in Sec. 3.1 of J. Santiago, J. Feng, S. Schuster, and M. Visser, *Immortality through the dark forces: Dark-charge primordial black holes as dark matter candidates* ([arXiv:2503.20696](https://arxiv.org/abs/2503.20696)). The scripts included in this repository may be regarded as supplementary material for the latter.

The code requires the 
following Julia packages:

    OrdinaryDiffEq, ArbNumerics, ImplicitEquations, LaTeXStrings, LsqFit, Plots, Pyplot
    
To run the code, open a terminal and cd to the following directory:

    > cd src/scripts/

To run the HWscript_DE.jl script, run the following command:

    > julia HWscript_DE.jl

Alternatively, you can simply paste the contents of HWscript_DE.jl into the Julia REPL (make sure Julia is running in the "src/scripts/" directory). This script will generate plots (describing the evolution of the mass and charge) which may be found in the "plots" directory.

Before running the script Contour_Plot.jl, please run Contour_Plot_Datagen.jl first (which supports multithreading).
