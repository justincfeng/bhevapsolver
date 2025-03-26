#-----------------------------------------------------------------------
#   Packages
#-----------------------------------------------------------------------

using OrdinaryDiffEq, Plots, ArbNumerics, LaTeXStrings

#-----------------------------------------------------------------------
#   Types
#-----------------------------------------------------------------------

setprecision(ArbFloat, 200)
tpfl    = ArbFloat{200}

#-----------------------------------------------------------------------
#   Condition for stopping integration
#-----------------------------------------------------------------------

function cond(u,t,integrator)   # Stop when Y hits 1e-6
    norm(u) - tpfl(1e-8)
end

effect!(integrator) = terminate!(integrator)

cb = ContinuousCallback(cond,effect!)

#-----------------------------------------------------------------------
#   Integrator
#-----------------------------------------------------------------------

integrator = AutoTsit5(Rosenbrock23())

rt  = tpfl("1e-16")     # Relative tolerance
at  = tpfl("1e-16")     # Absolute tolerance

