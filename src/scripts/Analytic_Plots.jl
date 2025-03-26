#-----------------------------------------------------------------------
#
#   ANALYTIC PLOTS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Packages
#-----------------------------------------------------------------------

using Plots
using ArbNumerics, LaTeXStrings, ImplicitEquations
using OrdinaryDiffEq

#-----------------------------------------------------------------------
#   Plot Settings
#-----------------------------------------------------------------------

default(
    size=(340, 240),
    #dpi=300,
    fontfamily="Times",
    tickfontsize=8,
    legendfontsize=9,
    guidefontsize=9,
    margin=3Plots.mm,
    legend=:topright
)

ylms=(0,1.05)
c1=:steelblue
c2=:red

labelphase = ( L"\mu" , L"Y" )

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
#   Type Definitions
#-----------------------------------------------------------------------

setprecision(ArbFloat, 200)
tpfl = ArbFloat{200}      

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
#
#   PLOTS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   ALLOWED REGION PHASE SPACE PLOT
#-----------------------------------------------------------------------

plot(range(0.9,3,length=1000),μ->HWslv.YH(μ,(0.9)^2), 
    fillrange = μ->1, 
    fillalpha = 0.20, 
    c = 1, lw=0,fillcolor=:grey50,linecolor=:grey90,
    label = L"t_\mathrm{evap} > t_\mathrm{univ}",
    xlabel=labelphase[1],ylabel=labelphase[2],
    fontfamily="Times",
    aspect_ratio = :equal, 
    yticks=([1.0], ["1"]), 
    xticks=([p[2],0.9], 
    [L"z_0",L"\mu_\mathrm{univ}"]), 
    size=(625,225),
    grid=false,framestyle=:axes,legend=:right, 
    margin=1Plots.mm)
plot!([-0.088, 3.01],[0, 0],arrow=(:closed, 0.1),color=:black, 
    label="",linewidth=1)  
plot!([0, 0], [0.2, 1.11], arrow=(:closed, 0.1), color=:black, label="",
    linewidth=1) 
plot!(μ->tpfl(1), 0, 3, xlims=(0, 3), ylims=(0, 1.1), label="", 
    linecolor=:grey70, linestyle=:dot)
plot!(Y->HWslv.muattr(Y,p),Y->Y,tpfl(0),tpfl(1),xlims=(0, 3),
    ylims=(0, 1.1),label=L"Y_H(\mu)",linecolor=c1)
plot!(s->0.9,s->s,tpfl(0.005),tpfl(0.995),xlims=(0, 3),ylims=(0, 1.1),
    label="",linestyle=:dash,linecolor=:orange)
plot!(s->p[2],s->s,tpfl(0.005),tpfl(0.997),xlims=(0, 3),ylims=(0, 1.1),
    label="",linestyle=:dashdot,linecolor=:grey50)
plot!(μ->HWslv.YH(μ,(0.9)^2),tpfl(0.9),tpfl(3),xlims=(0, 3),
    ylims=(0, 1.1),label="",linecolor=c1)
plot!(μ->tpfl(1),p[2],tpfl(0.9),xlims=(0, 3),ylims=(0, 1.1),
    label="",linecolor=c1)
savefig(dir*"PhaseSpace_AllowedRegion.pdf")

#-----------------------------------------------------------------------
#   PHASE SPACE mubar
#-----------------------------------------------------------------------

plot(range(0.9,3,length=1000),μ->HWslv.YH(μ,(0.9)^2), 
    fillrange = μ->1, 
    fillalpha = 0.20, 
    c = 1, lw=0,fillcolor=:grey50,linecolor=:grey90,
    label = L"t_\mathrm{evap} > t_\mathrm{univ}",
    xlabel=labelphase[1],ylabel=labelphase[2],
    fontfamily="Times",
    aspect_ratio = :equal, 
    yticks=([1.0], ["1"]), 
    xticks=([0.62*p[2],p[2],0.9], 
    [L"\bar\mu_1",L"\bar\mu_2",L"\bar\mu_3"]), 
    size=(625,225),
    grid=false,legend=:right, 
    margin=1Plots.mm)
plot!(range(0.9,3,length=1000),μ->HWslv.YH(μ,(0.9)^2), 
    lw=1,linecolor=:skyblue2,
    label = L"\bar\mu_3>z_0")
plot!(s->p[2],s->s,tpfl(0.005),tpfl(0.997),xlims=(0, 3),ylims=(0, 1.1),
    label="",linestyle=:dot,linecolor=:firebrick1)
plot!(μ->tpfl(1), 0, 3, xlims=(0, 3), ylims=(0, 1.1), label="", 
    linecolor=:grey70, linestyle=:dot)
plot!(s->0.9,s->s,tpfl(0.005),tpfl(0.995),xlims=(0, 3),ylims=(0, 1.1),
    label="",linestyle=:dot,linecolor=:skyblue2)
plot!(μ->HWslv.YH(μ,(0.9)^2),tpfl(0.9),tpfl(3),xlims=(0, 3),
    ylims=(0, 1.1),label="",linecolor=:skyblue2)
plot!(μ->tpfl(1),p[2],tpfl(0.9),xlims=(0, 3),ylims=(0, 1.1),
    label="",linecolor=:skyblue2)
plot!(μ->HWslv.YH(μ,(p[2])^2),tpfl(p[2]),tpfl(3),lw=1,
    linecolor=:firebrick1,label = L"\bar\mu_2=z_0")
plot!(Y->HWslv.muattr(Y,p),Y->Y,tpfl(0),tpfl(1),xlims=(0, 3),lw=1,
    label="",ylims=(0, 1.1),linecolor=:firebrick1)
plot!(s->tpfl(0.62*p[2]),s->s,tpfl(0.005),tpfl(0.935),xlims=(0, 3),
    ylims=(0, 1.1),label="",linestyle=:dot,linecolor=:black)
plot!(μ->HWslv.YH(μ,(6*p[2]/10)^2),tpfl(0.62*p[2]),tpfl(2.95),lw=1,
    linecolor=:black,label = L"\bar\mu_1<z_0")
plot!(Y->HWslv.muattr(Y,p),Y->Y,tpfl(0),tpfl(0.935),xlims=(0, 3),lw=1,
    label="",linecolor=:black,framestyle=:axes)
plot!([-0.088, 3.01], [0, 0], arrow=(:closed, 0.1), color=:black, 
    label="", linewidth=1)  
plot!([0, 0], [0.2, 1.11], arrow=(:closed, 0.1), color=:black, 
    label="", linewidth=1)
plot!([0.9, 3], [1, 1], linestyle=:dot, linecolor=:magenta, 
    label="", linewidth=2)

savefig(dir*"PhaseSpace_mubarS.pdf")

#-----------------------------------------------------------------------
#   ANALYTIC CONFIGURATION SPACE SOLUTION
#-----------------------------------------------------------------------

plot(μ->HWslv.YH(μ,(0.7)^2),
    tpfl(0.7),tpfl(1),
    xlims=(0.7, 1),ylims=ylms,
    label=L"Y_H(\mu)",xlabel=labelphase[1],ylabel=labelphase[2],
    fontfamily="Times",framestyle=:box,linecolor=c1,
    legend=(0.25,0.30))
plot!(μ->tpfl(1),p[2],tpfl(0.7),xlims=(0, 1),ylims=ylms,
    label="",linecolor=c1)
plot!(Y->HWslv.muattr(Y,p),Y->Y,tpfl(0),tpfl(0.98),xlims=(0, 1),
    ylims=ylms,label="",linecolor=c1)
plot!(Y->HWslv.muattr(Y,p),Y->Y,tpfl(0.98),tpfl(1),xlims=(0, 1),
    ylims=ylms,label="",linecolor=c1)
plot!(s->0.7,s->s,tpfl(0),tpfl(0.995),xlims=(0, 1),ylims=ylms,
    label=L"\mu_h",linestyle=:dash,linecolor=:orange)
plot!(s->p[2],s->s,tpfl(0),tpfl(0.997),xlims=(0, 1),ylims=ylms,
    label=L"z_0",linestyle=:dot,linecolor=:grey50)

savefig(dir*"PhaseSpace_Analytic.pdf")

#-----------------------------------------------------------------------
#   CONFIGURATION SPACE CURVE PLOT PARAMETERS
#-----------------------------------------------------------------------

plot_margin = 3Plots.mm
plot_size   = (340, 240)

Y1v = [ tpfl(0.003) , tpfl(0.02) , tpfl(0.07) , tpfl(0.18)  ]          
μcH = [ tpfl(0.0940891) , tpfl(0.173539) , tpfl(0.280817) , 
        tpfl(0.426952)  ]                                              
μ0v = [ tpfl(0.1) , tpfl(0.2) , tpfl(0.3) , tpfl(0.4)  ]               
μcS = [ tpfl(0.0517711) , tpfl(0.114169) , tpfl(0.197821) , 
        tpfl(0.313353)  ]                                              

#-----------------------------------------------------------------------
#   CONFIGURATION SPACE CURVE - MASS DISSIPATION ZONE
#-----------------------------------------------------------------------

plot(Y->HWslv.muattr(Y,p),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    legend=false,
    xlabel=labelphase[1],
    ylabel=labelphase[2],
    fontfamily="Times",
    framestyle=:box,
    linecolor=c2,
    size=plot_size,
    margin=plot_margin)
plot!(μ->HWslv.YH(μ,Y1v[1]),μcH[1],tpfl(1),xlims=(0, 1),ylims=ylms,
    linecolor=c1)
plot!(μ->HWslv.YH(μ,Y1v[1]),√(Y1v[1]),μcH[1],xlims=(0, 1),ylims=ylms,
    linestyle=:dash,linecolor=c1)
plot!(μ->HWslv.YH(μ,Y1v[2]),μcH[2],tpfl(1),xlims=(0, 1),ylims=ylms,
    linecolor=c1)
plot!(μ->HWslv.YH(μ,Y1v[2]),√(Y1v[2]),μcH[2],xlims=(0, 1),ylims=ylms,
    linestyle=:dash,
    linecolor=c1)
plot!(μ->HWslv.YH(μ,Y1v[3]),μcH[3],tpfl(1),xlims=(0, 1),ylims=ylms,
    linecolor=c1)
plot!(μ->HWslv.YH(μ,Y1v[3]),√(Y1v[3]),μcH[3],xlims=(0, 1),ylims=ylms,
    linestyle=:dash,
    linecolor=c1)
plot!(μ->HWslv.YH(μ,Y1v[4]),μcH[4],tpfl(1),xlims=(0, 1),ylims=ylms,
    linecolor=c1)
plot!(μ->HWslv.YH(μ,Y1v[4]),√(Y1v[4]),μcH[4],xlims=(0, 1),ylims=ylms,
    linestyle=:dash,
    linecolor=c1)
plot!(μ->HWslv.YH(μ,0.6^2),0.6,tpfl(1),xlims=(0, 1),ylims=ylms,
    linecolor=c1)
plot!(μ->HWslv.YH(μ,0.8^2),0.8,tpfl(1),xlims=(0, 1),ylims=ylms,
    linecolor=c1)

savefig(dir*"PhaseSpace_HawkingAttr.pdf")

#-----------------------------------------------------------------------
#   CONFIGURATION SPACE CURVE - CHARGE DISSIPATION ZONE
#-----------------------------------------------------------------------

plot(Y->HWslv.muattr(Y,p),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    legend=false,
    xlabel=labelphase[1],
    ylabel=labelphase[2],
    fontfamily="Times",
    framestyle=:box,
    linecolor=c2,
    size=plot_size,
    margin=plot_margin)
plot!(μ->HWslv.YS(μ,μ0v[1]),μcS[1],μ0v[1],xlims=(0, 1),ylims=ylms,
    linecolor=c1)
plot!(μ->HWslv.YS(μ,μ0v[1]),tpfl(0),μcS[1],xlims=(0, 1),ylims=ylms,
    linestyle=:dash,
    linecolor=c1)
plot!(μ->HWslv.YS(μ,μ0v[2]),μcS[2],μ0v[2],xlims=(0, 1),ylims=ylms,
    linecolor=c1)
plot!(μ->HWslv.YS(μ,μ0v[2]),tpfl(0),μcS[2],xlims=(0, 1),ylims=ylms,
    linestyle=:dash,
    linecolor=c1)
plot!(μ->HWslv.YS(μ,μ0v[3]),μcS[3],μ0v[3],xlims=(0, 1),ylims=ylms,
    linecolor=c1)
plot!(μ->HWslv.YS(μ,μ0v[3]),tpfl(0),μcS[3],xlims=(0, 1),ylims=ylms,
    linestyle=:dash,
    linecolor=c1)
plot!(μ->HWslv.YS(μ,μ0v[4]),μcS[4],μ0v[4],xlims=(0, 1),ylims=ylms,
    linecolor=c1)
plot!(μ->HWslv.YS(μ,μ0v[4]),tpfl(0),μcS[4],xlims=(0, 1),ylims=ylms,
    linestyle=:dash,
linecolor=c1)

savefig(dir*"PhaseSpace_SchwingerAttr.pdf")

#-----------------------------------------------------------------------
#   SHIFT IN ATTRACTOR DUE TO CHANGES IN MASS
#-----------------------------------------------------------------------

plot(Y->HWslv.muattr(Y,RSHWgen(0.8,1.0,1.0,tpfl)),Y->Y,
    tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    xlabel=labelphase[1],ylabel=labelphase[2],
    fontfamily="Times",framestyle=:box,
    linecolor=:lightsteelblue2,
    label=L"m_\chi=0.8m_{e}",
    legend=:bottomright)
plot!(Y->HWslv.muattr(Y,RSHWgen(0.9,1.0,1.0,tpfl)),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    linecolor=:lightsteelblue3,
    label=L"m_\chi=0.9m_{e}")
plot!(Y->HWslv.muattr(Y,RSHWgen(1.0,1.0,1.0,tpfl)),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    linecolor=c2,
    label=L"m_\chi=m_{e}")
plot!(Y->HWslv.muattr(Y,RSHWgen(1.1,1.0,1.0,tpfl)),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    linecolor=:steelblue3,
    label=L"m_\chi=1.1m_{e}")
plot!(Y->HWslv.muattr(Y,RSHWgen(1.2,1.0,1.0,tpfl)),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    linecolor=:steelblue4,
    label=L"m_\chi=1.2m_{e}")

savefig(dir*"mass_attractor.pdf")

#-----------------------------------------------------------------------
#   SHIFT IN ATTRACTOR DUE TO CHANGES IN CHARGE
#-----------------------------------------------------------------------

plot(Y->HWslv.muattr(Y,RSHWgen(1.0,0.8,1.0,tpfl)),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    xlabel=labelphase[1],ylabel=labelphase[2],
    fontfamily="Times",framestyle=:box,
    linecolor=:lightsteelblue2,
    label=L"e_\chi=0.8e",
    legend=:bottomright)
plot!(Y->HWslv.muattr(Y,RSHWgen(1.0,0.9,1.0,tpfl)),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    linecolor=:lightsteelblue3,
    label=L"e_\chi=0.9e")
plot!(Y->HWslv.muattr(Y,RSHWgen(1.0,1.0,1.0,tpfl)),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    linecolor=c2,
    label=L"e_\chi=e")
plot!(Y->HWslv.muattr(Y,RSHWgen(1.0,1.1,1.0,tpfl)),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    linecolor=:steelblue3,
    label=L"e_\chi=1.1e")
plot!(Y->HWslv.muattr(Y,RSHWgen(1.0,1.2,1.0,tpfl)),Y->Y,tpfl(0),tpfl(1),
    xlims=(0, 1),ylims=ylms,
    linecolor=:steelblue4,
    label=L"e_\chi=1.2e")

savefig(dir*"charge_attractor.pdf")

#-----------------------------------------------------------------------
#   RESCALED TIME INTEGRAL PLOT
#-----------------------------------------------------------------------

plot(μ->Δτatt(μ,z0),zero(tpfl),z0,
    legend=(0.15,0.6),
    xlabel=labelphase[1],
    ylabel=L"\Delta\tau_{\mathrm{att}}",
    label=L"\mathrm{Integrated}",
    fontfamily="Times",
    framestyle=:box,
    linecolor=c1,
    size=plot_size,
    margin=plot_margin)
plot!(τfit,zero(tpfl),z0,
    label=L"\mathrm{Fitted}",
    fontfamily="Times",
    linestyle=:dot,
    linecolor=:blue)
hline!([Δτatt(z0,z0)], 
    linestyle=:dash, 
    linecolor=c2, 
    alpha=0.5,
    label=L"0.30547"
    )

savefig(dir*"dtau_integral.pdf")

plot(μ->abs((τfit(μ)-Δτatt(μ,z0))/Δτatt(μ,z0)),tpfl(1e-4),z0,
    legend=false,
    xscale=:log10,
    yscale=:log10,
    xticks = 10.0 .^ (-4:1:0),
    yticks = 10.0 .^ (-6:-2),
    xlabel=labelphase[1],
    ylabel=L"|\Delta\tau_{\mathrm{att,f}}-\Delta\tau_{\mathrm{att,n}}|/\Delta\tau_{\mathrm{att,n}}",
    fontfamily="Times",
    framestyle=:box,
    linecolor=c1,
    size=plot_size,
    margin=plot_margin)

savefig(dir*"dtau_fit_residual.pdf")

#-----------------------------------------------------------------------
#   
#-----------------------------------------------------------------------

μh = [ tpfl(0.70) , tpfl(0.75) , tpfl(0.80) , tpfl(0.85) , tpfl(0.90) ];

ϑlist   = ( tpfl(3.712697401710582e-43) , tpfl(3.25978667490128e-46)  ,
            tpfl(2.805592079419033e-49) , tpfl(2.3496402985065857e-52),
            tpfl(1.8911314550422192e-55)    )                          ;

σMlist  = ( tpfl(3.712697401710582e-25) , tpfl(3.25978667490128e-27)  ,
            tpfl(2.805592079419033e-29) , tpfl(2.349640298506586e-31) ,
            tpfl(1.8911314550422194e-33)    )                          ;

tevap = tpfl[7.565597362277192e31  4.0165595291330044e36  2.133024747097304e41  1.1327693380596368e46  6.0157141498241565e50    ;
6.45150462056052e23	9.075826324379927e27	1.2782284810050254e32	1.8003375203141252e36	2.5357150174058666e40;
5.233187493772178e15	1.9387217678221095e19	7.2125107374083895e22	2.684044500472054e26	9.98855051474276e29;
3.9891680558945365e7	3.8390725120363434e10	3.748554560146054e13	3.666209827963917e16	3.586338930250767e19;
0.28357224771886547	68.49678357177582	17309.943871379353	4.4102057753856005e6	1.1252447007940865e9] ;

TeLg = log10.(tevap)

jloop = 1:5

#-----------------------------------------------------------------------
for j=jloop
    ix  = j+13
    ϑ   = ϑlist[j]
    σM  = σMlist[j]

    ξ   = (b0*z0*ϑ/(b0*z0+TWO*log(ϑ)))^(1/3)
    σm  = ϑ/(σM*ξ^2)
    σet = σm*(ϑ^(1/3))

    Ms0 = σM*(1e8)
    Ms  = tpfl(1474.07)*Ms0

    plot(   μ->log10(Δts(μ,σm,σet,tpfl)) , 0.65 , 0.95 , 
            xlabel= snsMh(Ms0),
            ylabel= L"\log_{10}(t_\mathrm{evap}/{\mathrm{yr}})",
            label="Estimate",
            fontfamily="Times",
            linecolor=:steelblue,
            framestyle=:box,
            legend=:bottomright,
            size=plot_size,
            margin=plot_margin,
            legendfontsize=6,
            markersize=3
        )
    scatter!(μh,log10.(tevap[j,:]),
            label="Numerical ",
            fontfamily="Times",
            markersize=3,
            markerstrokewidth=0,
            markeralpha=1)
    savefig(dir*"TimeVsM0_"*string(ix)*".pdf")
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------