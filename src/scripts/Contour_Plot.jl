#-----------------------------------------------------------------------
#
#   CONTOUR PLOTS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Packages
#-----------------------------------------------------------------------

using Plots
import PyPlot
using ArbNumerics, LaTeXStrings
using OrdinaryDiffEq
using Serialization

#-----------------------------------------------------------------------
#   Plot Settings
#-----------------------------------------------------------------------

PyPlot.rc("text", usetex=true)
PyPlot.rc("text.latex", preamble=raw"\usepackage{amsmath}")
pyplot()

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
    ϑ = (σet/σm)^3
    ξ = (b0*z0*ϑ/(b0*z0+TWO*log(ϑ)))^(1/3)
    σM = ϑ/(σm*ξ^2)
    MS = σM * Ms / Msol
    τu = τuniv(σm,σet,tpfl)
    return MS*maximum([μfit(τu,tpfl),μfromτ(τu,B0,z0,Δτattz0,tpfl)])
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    curvMC1( σe , tpfl=Float64 )
    
This function computes the curve defined by the constraint 'MC3=Munv', 
assuming that 'Munv' is at its asymptotc value for small 'μ'.
"""
function curvMC1( σe , tpfl=Float64 )
    k   = tpfl(1e-15)
    ONE = one(tpfl)
    return k/Munv(ONE,ONE)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    MzMunivCrv(ϑ,tpfl=Float64,αn=-1)

This function computes the Mz=Muniv curve parametrically, as a function 
of the parameter ϑ.
"""
function MzMunivCrv(ϑ,Δτattz0=tpfl(0.3054715988206478),tpfl=Float64,
                    αn=-1)
    p = pars(tpfl)
    (s0,z0,b0)      = par[1]
    (qe,me,ħ)       = par[2]
    (Ms,Msol,tuniv) = par[3]
    (α1,α2)         = par[4]
    (c,yr2s,yr2m)   = par[6]
    (ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI) = par[5]

    α = α1

    tunv = tuniv

    ξ = (b0*z0*ϑ/(b0*z0+TWO*log(ϑ)))^(1/3)

    return (ϑ/ξ^2)*((Ms*Δτattz0/(s0*tunv))^(1/3)) .* [ONE,ONE*ξ]
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    tunivcurve(ϑ,tpfl=Float64,Δτatt=tpfl(0.3054715988206478))

This function computes the tuniv curve parametrically as a function of
the parameter ϑ.
"""
function tunivcurve(ϑ,tpfl=Float64,Δτatt=tpfl(0.3054715988206478))
    p = pars(tpfl)
    (s0,z0,b0)      = par[1]
    (qe,me,ħ)       = par[2]
    (Ms,Msol,tuniv) = par[3]
    (α1,α2)         = par[4]
    (c,yr2s,yr2m)   = par[6]
    ZER             = zero(tpfl)
    (ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI) = par[5]
    THRT = tpfl(30)

    α = α1

    tunv = tuniv

    ξ = (b0*z0*ϑ/(b0*z0+TWO*log(ϑ)))^(1/3)

    σm = THRT*Ms*((THRT*PI*Δτatt*ϑ*(b0*z0+TWO*log(ϑ))^2) / 
            (b0^2*tunv*z0^2*α*ħ))^(1/3)

    return σm .* [ONE,ONE*ξ]
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   PLOTS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   MINIMUM MASS CONTOUR PLOT
#-----------------------------------------------------------------------

(Mmn,τz,μz,C2,C3) = deserialize(dir * "data.jls")

x   = tpfl.( 10 .^ range(0, 12, length=4000) )
y   = tpfl.( 10 .^ range(-9, 0, length=2200) )

tv_all = collect(-25:1:-15)
tv = tv_all[2:1:end]
tl = [L"10^{%$i}" for i in tv]

p = Plots.contourf(log10.(x), log10.(y), log10.(Mmn),
    xlims=(0, 12), ylims=(-8, 0), levels=tv_all,
    colorbar=true, colorbar_ticks=(tv, tl),
    color=cgrad(:Spectral_11,scale = :log),
    background_color_inside=:black, legend=(0.34,0.11),
    grid=false, framestyle=:box, aspect_ratio=:equal, 
    title=L"M_\mathrm{min}/M_{\odot}", 
    xlabel=L"\log_{10}(m_\chi/m_e)", 
    ylabel=L"\log_{10}(e_\chi/e)",
    tickfontsize=16, guidefontsize=18, titlefontsize=20,
    colorbar_tickfontsize=16, legendfontsize=11)

plot!(x->log10.(MzMunivCrv(x)[1]),x->log10.(MzMunivCrv(x)[2]),
    tpfl(1e-55),tpfl(1e-40),label=L"M_{\mathrm{univ}}=M_{z_0}",
    linewidth=1, color=:black)

plot!(x->log10.(curvMC1(x)),x->log10.(x),tpfl(1e-8),
    tpfl(1.0),label=L"M_{\mathrm{univ}}\approx\hbar/m_\mathrm{e}",
    linewidth=1, color=:black,linestyle=:dot)

Plots.savefig(dir * "contour_plot_Muniv.pdf")

#-----------------------------------------------------------------------
#   ATTRACTOR TIME PLOT
#-----------------------------------------------------------------------

tv_all = collect(-10:5:95)
tv = tv_all[1:2:end]
tl = [L"10^{%$i}" for i in tv]

p = Plots.contourf(log10.(x), log10.(y), log10.(τz),
    xlims=(0, 12), ylims=(-8, 0), levels=tv_all,
    colorbar=true, colorbar_ticks=(tv, tl), color=:Spectral_11,
    background_color_inside=:black, legend=:topleft,
    grid=false, framestyle=:box, aspect_ratio=:equal, 
    title=L"\Delta t_\mathrm{att}(M_{z_0})/\mathrm{yr}", 
    xlabel=L"\log_{10}(m_\chi/m_e)", 
    ylabel=L"\log_{10}(e_\chi/e)",
    tickfontsize=16, guidefontsize=18, titlefontsize=20,
    colorbar_tickfontsize=16, legendfontsize=11)

plot!(x->log10.(tunivcurve(x)[1]),x->log10.(tunivcurve(x)[2]),
    tpfl(1e-65),tpfl(1e-42),
    label=L"\Delta t_\mathrm{att}=t_\mathrm{univ}",
    linewidth=1, color=:black,linestyle=:dash)

plot!(x->log10.(MzMunivCrv(x)[1]),x->log10.(MzMunivCrv(x)[2]),
    tpfl(1e-54),tpfl(1e-40),label=L"M_{\mathrm{univ}}=M_{z_0}",
    linewidth=1, color=:black)

Plots.savefig(dir * "contour_plot_delta_t_att.pdf")

#-----------------------------------------------------------------------
#   RESCALED MASS PLOT
#-----------------------------------------------------------------------

tv_all = collect(-27:1:1)
tv = tv_all[1:3:end]
tl = [L"10^{%$i}" for i in tv]

p = Plots.contourf(log10.(x), log10.(y), log10.(μz),
    xlims=(0, 12), ylims=(-8, 0), levels=tv_all,
    colorbar=true, colorbar_ticks=(tv, tl), color=:Spectral_11,
    background_color_inside=:black, legend=:bottomleft,
    grid=false, framestyle=:box, aspect_ratio=:equal, 
    title=L"\mu_\mathrm{univ}/z_0", 
    xlabel=L"\log_{10}(m_\chi/m_e)", 
    ylabel=L"\log_{10}(e_\chi/e)",
    tickfontsize=16, guidefontsize=18, titlefontsize=20,
    colorbar_tickfontsize=16, legendfontsize=11)

plot!(x->log10.(MzMunivCrv(x)[1]),x->log10.(MzMunivCrv(x)[2]),
    tpfl(1e-54),tpfl(1e-40),label=L"M_{\mathrm{univ}}=M_{z_0}",
    linewidth=1, color=:black)

Plots.savefig(dir * "contour_plot_muz.pdf")

#-----------------------------------------------------------------------
