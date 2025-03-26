using DelimitedFiles
using LsqFit
using Plots, LaTeXStrings
using OrdinaryDiffEq

#-----------------------------------------------------------------------
#   Parameters
#-----------------------------------------------------------------------

tpfl    = Float64

#-----------------------------------------------------------------------
#   Includes
#-----------------------------------------------------------------------

include("../../HWslv.jl")
include("../../HWFunc.jl")

dir     ="../../../plots/"
if !isdir(dir)
    mkpath(dir)
end

#-----------------------------------------------------------------------
#   Parameters
#-----------------------------------------------------------------------

p = pars(tpfl)

(s0,z0,b0)          = p[1] 
(qe,me,ħ)           = p[2]
(Ms,Msol,tuniv)     = p[3]
(α1,α2)             = p[4]
(c,yr2s,yr2m)       = p[6]
(ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI) = p[5]

#-----------------------------------------------------------------------
#   Model
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    model(x, p)

The model we use for the fit is a simplified rational function of the
form: 'x^p1/(p2*x^p1 + p3)'.
"""
function model(x, p)
    (x .^ p[1]) ./ (p[2] .* x .^ p[1] .+ p[3])
end  #------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    imodel(x, p)

The inverse of the model used for the fit has the following form: 
'(p3*x/(1-p2*x^p1))^(1/p1)'.
"""
function imodel(x, p)
    ((x .* p[3]) ./ (one(typeof(x)) .- p[2] .* x))^(1/p[1])
end  #------------------------------------------------------------------

# Reduced function
Δτatt0 = μ->Δτatt(μ,z0)

# Domain breaks
μ1 = 0.07
μ2 = 0.22
μ3 = 0.35

τ1 = Δτatt0(μ1)
τ2 = Δτatt0(μ2)
τ3 = Δτatt0(μ3)

# Evaluated points
npoints = 30000

μvals1 = vcat(  tpfl.(range(μ1, μ1*1.01, length=npoints)),
                tpfl.(range(μ1, μ2, length=npoints))    ,
                tpfl.(range(μ2, μ2*0.99, length=npoints)))
τvals1 = Δτatt0.(μvals1)

μvals2 = vcat(  tpfl.(range(μ2, μ2*1.01, length=npoints)),
                tpfl.(range(μ2, μ3, length=npoints))    ,
                tpfl.(range(μ3, μ3*0.99, length=npoints)))
τvals2 = Δτatt0.(μvals2)

μvals3 = vcat(  tpfl.(range(μ3, μ3*1.01, length=npoints)),
                tpfl.(range(μ3, z0, length=npoints))    ,
                tpfl.(range(z0, z0*0.99, length=npoints)))
τvals3 = Δτatt0.(μvals3)

#-----------------------------------------------------------------------
#   Fit
#-----------------------------------------------------------------------

pi1  = [5.0, 1.0, 1.0]
pi2  = [5.0, 1.0, 1.0]
pi3  = [5.0, 1.0, 1.0]

fit1 = curve_fit(model, μvals1, τvals1, pi1)
fit2 = curve_fit(model, μvals2, τvals2, pi2)
fit3 = curve_fit(model, μvals3, τvals3, pi3)

#-----------------------------------------------------------------------
"""
    Δτatt1(x,p1=fit1.param,p2=fit2.param,p3=fit3.param,x1=μ1,x2=μ2,
             x3=μ3,x4=z0)

Here is the fitted function.
"""
function Δτatt1(x,p1=fit1.param,p2=fit2.param,p3=fit3.param,x1=μ1,x2=μ2,
                x3=μ3,x4=z0)
    tpfl=typeof(x)
    if x < x1
        return 32*x^3/81
    elseif x < x2
        return model(x,p1)
    elseif x <= x3
        return model(x,p2)
    elseif x <= x4
        return model(x,p3)
    elseif x >= x4
        return model(x4,p3)
    else
        return zero(tpfl)
    end
end     #---------------------------------------------------------------

discont1 =  100*(Δτatt1(μ1*(1.0001))-Δτatt1(μ1*(0.9999))) / 
                (Δτatt1(μ1*(1.0001))+Δτatt1(μ1*(0.9999)))

discont2 =  100*(Δτatt1(μ2*(1.0001))-Δτatt1(μ2*(0.9999))) / 
                (Δτatt1(μ2*(1.0001))+Δτatt1(μ2*(0.9999)))

discont3 =  100*(Δτatt1(μ3*(1.0001))-Δτatt1(μ3*(0.9999))) / 
                (Δτatt1(μ3*(1.0001))+Δτatt1(μ3*(0.9999)))

open(dir * "fit_results.txt", "w") do io
    println(io, "# Model: x^p1/(p2*x^p1 + p3)\n")
    println(io, "# Model Fit")
    println(io, "# μ ∈ ($μ1,$μ2)")
    println(io, "Params: ", fit1.param)
    println(io, "Std Errors: ", sum(fit1.resid)/npoints)
    println(io, "R-squared: ", 1 - sum(fit1.resid.^2)/
                               sum((τvals1 .- sum(τvals1)/npoints).^2))
    println(io, "Discontinuity at μ = $μ1: ", abs(discont1), "%")
    
    println(io, "\n# μ ∈ ($μ2,$μ3)")
    println(io, "Params: ", fit2.param)
    println(io, "Std Errors: ", sum(fit2.resid)/npoints)
    println(io, "R-squared: ", 1 - sum(fit2.resid.^2)/
                               sum((τvals2 .- sum(τvals2)/npoints).^2))
    println(io, "Discontinuity at μ = $μ2: ", abs(discont2), "%")

    println(io, "\n# μ ∈ ($μ3,$z0)")
    println(io, "Params: ", fit3.param)
    println(io, "Std Errors: ", sum(fit3.resid)/npoints)
    println(io, "R-squared: ", 1 - sum(fit3.resid.^2)/
                               sum((τvals3 .- sum(τvals3)/npoints).^2))
    println(io, "Discontinuity at μ = $μ3: ", abs(discont3), "%")
end

#-----------------------------------------------------------------------
#   Plot Fit
#-----------------------------------------------------------------------

plot_margin = 3Plots.mm
plot_size = (340, 240)

plot(Δτatt1,0,z0,xlims=(0,1.1*z0),ylims=(0,0.32),xlabel=L"\mu",
    ylabel=L"Y",fontfamily="Times",linecolor=:steelblue,framestyle=:box,
    label="Fitted Function",legend=:topleft,size=plot_size,
    margin=plot_margin)

savefig(dir * "fitted_function.pdf")
