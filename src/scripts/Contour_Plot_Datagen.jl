#-----------------------------------------------------------------------
#
#   CONTOUR PLOT DATA GENERATION
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Packages
#-----------------------------------------------------------------------

using Plots
import PyPlot
using ArbNumerics, LaTeXStrings
using OrdinaryDiffEq
using Base.Threads
using Serialization

#-----------------------------------------------------------------------
#   Parameters
#-----------------------------------------------------------------------

tpfl    = Float64

#-----------------------------------------------------------------------
#   Includes
#-----------------------------------------------------------------------

include("../HWslv.jl")
include("../HWFunc.jl")
include("PlotFunc.jl")
include("Fitted.jl")

dir     ="../../plots/"
if !isdir(dir)
    mkpath(dir)
end

#-----------------------------------------------------------------------
#   PARAMETERS
#-----------------------------------------------------------------------

par = pars(tpfl)

(s0,z0,b0)      = par[1]
(qe,me,ħ)       = par[2]
(Ms,Msol,tuniv) = par[3]
(α1,α2)         = par[4]
(c,yr2s,yr2m)   = par[6]
ZER             = zero(tpfl)
(ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI) = par[5]

p   = pars(tpfl)[1]

#-----------------------------------------------------------------------
#   MASS FUNCTION
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    Munv( σm , σet , tpfl=Float64 )
    
This function computes the mass of a black hole with an evaporation time 
at the age of the universe if it is initially near extremal or
on the attractor curve.
"""
function Munv(σm,σet,tpfl=Float64)
    Δτattz0=tpfl(0.3054715988206478)
    B0 = b0+SIX*log(σet/σm)/z0
    #σM = (σet/σm^2)*(ONE+SIX*log(σet/σm)/(b0*z0))
    ϑ = (σet/σm)^3
    ξ = (b0*z0*ϑ/(b0*z0+TWO*log(ϑ)))^(1/3)
    σM = ϑ/(σm*ξ^2)
    MS = σM * Ms / Msol
    τu = τuniv(σm,σet,tpfl)
    if τu>Δτattz0
        return MS*μfromτ(τu,B0,z0,Δτattz0,tpfl)
    elseif τu<Δτattz0
        return MS*μfit(τu,tpfl)
    else
        return MS*z0
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    μunvz0( σm , σet , tpfl=Float64 )
    
This function computes the mass of a black hole with an evaporation time 
at the age of the universe if it is initially near extremal or
on the attractor curve.
"""
function μunvz0(σm,σet,tpfl=Float64)
    Δτattz0=tpfl(0.3054715988206478)
    B0 = b0+SIX*log(σet/σm)/z0
    ϑ = (σet/σm)^3
    ξ = (b0*z0*ϑ/(b0*z0+TWO*log(ϑ)))^(1/3)
    σM = ϑ/(σm*ξ^2)
    MS = σM * Ms / Msol
    τu = τuniv(σm,σet,tpfl)
    if τu>Δτattz0
        return μfromτ(τu,B0,z0,Δτattz0,tpfl)/z0
    elseif τu<Δτattz0
        return μfit(τu,tpfl)/z0
    else
        return ONE
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    Mmin( σm , σet , tpfl=Float64 )
    
This function computes the maximum of 'Munv' and 'MC1'.
"""
function Mmin(σm,σet,tpfl=Float64)
    return maximum([Munv(σm,σet,tpfl),MC1(σm,σet,tpfl)])
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    constraint2( σm , σet , tpfl=Float64 )
    
This function computes constraint ii.
"""
function constraint2( σm , σet , tpfl=Float64 )
    k   = qe^3/(Ms*z0*me^2)
    ONE = one(tpfl)
    f   = ONE - σet^2 * (b0*z0/(b0*z0+SIX*log(σet/σm)))^(5/3)
    return hvs(f)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    constraint3( σm , σet , tpfl=Float64 )
    
This function computes constraint iii.
"""
function constraint3( σm , σet , tpfl=Float64 )
    f   = b0*z0+TWO*log((σet/σm)^3)-ONE
    return hvs(f)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   PLOTS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   MINIMUM MASS CONTOUR PLOT
#-----------------------------------------------------------------------

x   = tpfl.( 10 .^ range(0, 12, length=4000) )
y   = tpfl.( 10 .^ range(-9, 0, length=2200) )
z   = zeros(tpfl, (length(y),length(x)) )

(Mmn,τz,μz,C2,C3) = (copy(z),copy(z),copy(z),copy(z),copy(z))

Threads.@threads for i=1:length(x)
    for j=1:length(y)
        C2[j,i]  = constraint2(x[i],y[j],tpfl)
        C3[j,i]  = constraint3(x[i],y[j],tpfl)
        Mmn[j,i] = C2[j,i] * C3[j,i] * Mmin(x[i],y[j],tpfl)
        τz[j,i]  = C2[j,i] * C3[j,i] * Δt(z0,x[i],y[j],tpfl)
        μz[j,i]  = C2[j,i] * C3[j,i] * μunvz0(x[i],y[j],tpfl)
    end
end

serialize(dir * "data.jls", (Mmn,τz,μz,C2,C3))

#-----------------------------------------------------------------------
