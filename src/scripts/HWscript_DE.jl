#-----------------------------------------------------------------------
#
#   HW SCRIPT with Dark Electomagnetism Parameters
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Includes
#-----------------------------------------------------------------------

include("../HWglobal.jl")
include("../HWslv.jl")
include("../HWFunc.jl")
include("PlotFunc.jl")

#-----------------------------------------------------------------------
#   Plot directories
#-----------------------------------------------------------------------

dir     ="../../plots/"
if !isdir(dir)
    mkpath(dir)
end

#-----------------------------------------------------------------------
#   Parameters
#-----------------------------------------------------------------------

par = pars(tpfl)

(s0,z0,b0)          = par[1]
(qe,me,ħ)           = par[2]
(Msc,Msol,tuniv)    = par[3]
(α1,α2)             = par[4]
(c,yr2s,yr2m)       = par[6]
ZER                 = zero(tpfl)
(ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI) = par[5]

μh = [ tpfl(0.70) , tpfl(0.75) , tpfl(0.80) , tpfl(0.85) , tpfl(0.90) ]

YIDlist = μh .^ 2 ;

ϑlist   = ( tpfl(3.712697401710582e-43) , tpfl(3.25978667490128e-46)  ,
            tpfl(2.805592079419033e-49) , tpfl(2.3496402985065857e-52),
            tpfl(1.8911314550422192e-55)    )

σMlist  = ( tpfl(3.712697401710582e-25) , tpfl(3.25978667490128e-27)  ,
            tpfl(2.805592079419033e-29) , tpfl(2.349640298506586e-31) ,
            tpfl(1.8911314550422194e-33)    )

α       = α1

pa = (s0,z0,b0) 

#-----------------------------------------------------------------------
#   Main loop
#-----------------------------------------------------------------------
 
tevap   = zeros(tpfl,(5,5))             ;
MsMat   = zeros(tpfl,(5,5))             ;
σmMat   = zeros(tpfl,(5,5))             ;
σeMat   = zeros(tpfl,(5,5))             ;
ξMat    = zeros(tpfl,(5,5))             ;
TFM     = zeros(tpfl,(5,5))             ;

jloop = 1:5                             ;
iloop = 1:5                             ;

for j=jloop, i=iloop

    ix  = j+13

    ϑ   = ϑlist[j]
    σM  = σMlist[j]

    ξ   = (b0*z0*ϑ/(b0*z0+TWO*log(ϑ)))^(1/3)
    σm  = ϑ/(σM*ξ^2)
    σet = σm*(ϑ^(1/3))

    σe  = ξ*σm

    (S0,Z0,B0,σMgen) = RSHW(σm,σet,tpfl)

    p = (S0,Z0,B0)

    Ms0 = σM*(1e8)
    Ms  = tpfl(1474.07)*Ms0

    MsMat[j,i] = Ms0
    σeMat[j,i] = σe
    σmMat[j,i] = σm
    ξMat[j,i]  = ξ

    u0  = [ tpfl(1) , YIDlist[i] ]

    f   = (u,p,t)->HWslv.F(u,p)

    tspan1 = tpfl.((0.0,2*Δτest(sqrt(YIDlist[i]),p)))  ;

    sol = solve(   ODEProblem(f,u0,tspan1,p),integrator,
                    reltol=rt,abstol=at,callback=cb)

    tf1  = sol.t[end]

    Xr1  = t->(sol(t)[1])
    Yr1  = t->(sol(t)[2])

    X1   = t->Xr1(t*tf1)
    Y1   = t->Yr1(t*tf1)

    X1s  = s->X1(ers(s))
    Y1s  = s->Y1(ers(s))

    ylms=(0,1.04)
    
    plot_margin = 3Plots.mm
    plot_size = (340, 240)
    
    c1=:steelblue3
    c2=:red3
    c3=:darkorange3

    labelphase = ( snsM(Ms0) , L"(Q/M)^2" )

    plot(s->μh[i],s->s,tpfl(0),tpfl(1),
        xlims=(0,1.02),ylims=ylms,
        size=plot_size,
        margin=plot_margin,
        label=L"\mu_h",
        linestyle=:dash,
        linecolor=c3,
        legend=(0.25,0.30),
        fontfamily="Times",
        framestyle=:box)
    plot!(s->p[2],s->s,tpfl(0),tpfl(1),xlims=(0,1.02),ylims=ylms,
        label=L"z_0",linestyle=:dot,linecolor=:grey30)
    plot!(X1s,Y1s,tpfl(0),tpfl(1e-2),xlims=(0,1.02),ylims=ylms,
        label=L"Y(\mu)",xlabel=labelphase[1],ylabel=labelphase[2],
        linecolor=c1)
    plot!(X1s,Y1s,tpfl(1e-2),tpfl(2e-2),xlims=(0,1.02),ylims=ylms,
        label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(2e-2),tpfl(3e-2),xlims=(0,1.02),ylims=ylms,
        label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(3e-2),tpfl(4e-2),xlims=(0,1.02),ylims=ylms,
        label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(4e-2),tpfl(6e-2),xlims=(0,1.02),ylims=ylms,
        label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(6e-2),tpfl(0.13),xlims=(0,1.02),ylims=ylms,
        label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(0.13),tpfl(0.87),xlims=(0,1.02),ylims=ylms,
        label="",linecolor=c2)
    plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
        ylims=ylms,label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
        ylims=ylms,label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
        ylims=ylms,label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
        ylims=ylms,label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),
        xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
        xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
        xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
        xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
    plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),
        ylims=ylms,label="",linecolor=c1)

    savefig(dir*"PhaseSpace_"*string(ix)*"_"*string(i)*".pdf")

    TF          = tf1
    TFyr        = Float64(TF*Ms/S0)*(1/(c*yr2s))
    TFM[j,i]    = TF

    tevap[j,i] = TFyr

    XFFn    = t->Xr1(TF*t)
    YFFn    = t->Yr1(TF*t)

    XFFns   = s->XFFn(ers(s))
    YFFns   = s->YFFn(ers(s))

    QFFn    = t->sqrt(YFFn(t))*XFFn(t)
    QFFns   = s->sqrt(YFFns(s))*XFFns(s)

    labelsolte = ( snst(TFyr) , snsM(Ms0) )

    plot(s->s,s->μh[i],tpfl(0),tpfl(1),
        xlims=(0, 1.02),ylims=ylms,
        size=plot_size,
        margin=plot_margin,
        label=L"\mu_h",
        linestyle=:dash,
        linecolor=c3,
        framestyle=:box)
    plot!(ers,XFFns,0,0.01,xlims=(0, 1.02),ylims=ylms,label="",
        xlabel=labelsolte[1],
        ylabel=labelsolte[2],
        fontfamily="Times",linecolor=c1,legend=:bottomleft)
    plot!(ers,XFFns,0.01,0.98,xlims=(0, 1.02),ylims=ylms,label="",
        fontfamily="Times",linecolor=c1)
    plot!(ers,XFFns,0.98,0.9999,xlims=(0, 1.02),ylims=ylms,
        fontfamily="Times",label=L"\mu(t)",linecolor=c1)

    savefig(dir*"TimeSoln_"*string(ix)*"_"*string(i)*".pdf")

    dYFFns = s->log10(one(tpfl)-√(YFFns(s)))

    plot(ers,dYFFns,0,0.01,xlims=(0, 1.02),label="",
        size=plot_size,
        margin=plot_margin,
        xlabel=labelsolte[1],
        ylabel=L"{\rm log}_{10}(1-\sqrt{Y})",
        fontfamily="Times",linecolor=c1)
    plot!(ers,dYFFns,0.01,0.98,xlims=(0, 1.02),label="",
        fontfamily="Times",linecolor=c1)
    plot!(ers,dYFFns,0.98,0.9999,xlims=(0, 1.02),
        fontfamily="Times",label="",linecolor=c1)

    savefig(dir*"dYSoln_"*string(ix)*"_"*string(i)*".pdf")

    print( string((j,i))*"\t"*string(Float64(μh[i]))*"\t" * 
           string(Float64(ξMat[j,i]))*"\t" * 
           string(Float64(sol.t[end]))*"\t"* 
           string(Float64(dYFFns(1/2)))*"\t" *
           string(Float64.(sol.u[end]))*"\n" )
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Parameters for time vs M0
#-----------------------------------------------------------------------
for j=jloop
    ix  = j+13
    ϑ   = ϑlist[j]
    σM  = σMlist[j]

    ξ   = (b0*z0*ϑ/(b0*z0+TWO*log(ϑ)))^(1/3)
    σm  = ϑ/(σM*ξ^2)
    σet = σm*(ϑ^(1/3))

    σe  = ξ*σm

    (S0,Z0,B0,σMgen) = RSHW(σm,σet,tpfl)

    p = (S0,Z0,B0)

    S0  = p[1]
    Ms0 = σM*(1e8)
    Ms  = tpfl(1474.07)*Ms0

    plot( x->log10(Δt(x,σm,σet,tpfl)) , tpfl(0.65) , tpfl(0.95), 
        xlabel= snsMh(Ms0),
        ylabel= L"\log_{10}(t_\mathrm{evap}/{\mathrm{yr}})",
        label="Estimate",fontfamily="Times",
        linecolor=:steelblue,framestyle = :box,
        legend=:bottomright
        )
    scatter!(μh,log10.(tevap[j,:]),label="Numerical",fontfamily="Times",
             markersize=2, legend_markersize=2, markerstrokewidth=0)
    savefig(dir*"TimeVsM0_"*string(ix)*".pdf")
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   WRITE OUTPUTS TO FILE
#-----------------------------------------------------------------------

(s0,z0,b0) = pa

file1   = open(dir*"DE_output.txt","w")
file2   = open(dir*"DE_table.txt","w")

write(file2, "EVAPORATION TIMES\n")

write(file2, "         𝜚/𝜚0         ")
for i=iloop
    write(file2, "\t         "*"μh="*string(Float64(μh[i]))*"         ")
end
write(file2, "\n")

for j=jloop
    write(file2, string(Float64(ξMat[j])))
    for i=iloop
        write(file2, "\t"*string(Float64(tevap[j,i])))
    end
    write(file2, "\n")
end

write(file2, "\n\n")

write(file2, "RELATIVE DIFFERENCES IN LOG OF EVAPORATION TIMES\n")

write(file2, "         𝜚/𝜚0         ")
for i=iloop
    write(file2, "\t         "*"μh="*string(Float64(μh[i]))*"         ")
end
write(file2, "\n")

for j=jloop
    write(file2, string(Float64(ξMat[j])))    
    for i=iloop
        ϑ   = ϑlist[j]
        σM  = σMlist[j]
        ξ   = (b0*z0*ϑ/(b0*z0+TWO*log(ϑ)))^(1/3)
        σm  = ϑ/(σM*ξ^2)
        σet = σm*(ϑ^(1/3))
        diff = abs(log10(Δts(μh[i],σm,σet,tpfl))-log10(tevap[j,i])) / 
                abs(log10(tevap[j,i]))
        write(file2, "\t"*string(Float64(diff)))
    end
    write(file2, "\n")
end

for j=jloop
    write(file1, "χ ="*"\t"*string(Float64.(ξMat[j,1]))*"\n")
    write(file1, "Ms ="*"\t"*string(Float64.(MsMat[j,:]))*"\n")
    write(file1, "σe ="*"\t"*string(Float64.(σeMat[j,:]))*"\n")
    write(file1, "σm ="*"\t"*string(Float64.(σeMat[j,:]))*"\n")
    write(file1, "\n")
    write(file1, "\n")
write(file1, "\n")
end

close(file1)
close(file2)
