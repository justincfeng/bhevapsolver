#-----------------------------------------------------------------------
#
#   HW SCRIPT (Hiscock and Weems's original parameters)
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
#   Plot directory
#-----------------------------------------------------------------------

dir     ="../../plots/"
if !isdir(dir)
    mkpath(dir)
end

#-----------------------------------------------------------------------
#   Parameters
#-----------------------------------------------------------------------

P   = pars(tpfl,0)

p0                  = P[1]
(qe,me,ħ)           = P[2]
(Ms,Msol,tuniv)     = P[3]
(α1,α2)             = P[4]
(c,yr2s,yr2m)       = P[6]
(ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI) = P[5]

α = α2

#-----------------------------------------------------------------------
#   FIGURE 3
#-----------------------------------------------------------------------

σM  = ONE

p   = p0

s0  = p[1]
Ms0 = σM*(1e8)
Ms  = tpfl(1474.07)*Ms0

u0  = [ tpfl(0.3) , tpfl(0.999) ]

f   = (u,p,t)->HWslv.F(u,p)

tspan1 = tpfl.((0.0,0.003))  ;

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

# Set consistent plot margins and spacing
plot_margin = 3Plots.mm
plot_size = (340, 240)

c1=:steelblue3
c2=:red3
c3=:darkorange3

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

plot(X1s,Y1s,tpfl(0),tpfl(1e-2),xlims=(0,1.02),ylims=ylms,label="",
    xlabel=labelphase[1],ylabel=labelphase[2],fontfamily="Times",
    linecolor=c1,framestyle=:box,size=plot_size,margin=plot_margin)
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)

    savefig(dir*"PhaseSpace_HW_3.pdf")

    TF          = tf1
    TFyr        = Float64(TF*Ms/s0)*(1/(3e8))*(1/3.154e7)

    XFFn    = t->Xr1(TF*t)
    YFFn    = t->Yr1(TF*t)

    XFFns   = s->XFFn(ers(s))
    YFFns   = s->YFFn(ers(s))

    QFFn    = t->sqrt(YFFn(t))*XFFn(t)
    QFFns   = s->sqrt(YFFns(s))*XFFns(s)

    labelsolte = ( snst(TFyr) , snsMu(Ms0,2) )

plot(ers,XFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 0.3),label="",
    xlabel=labelsolte[1],ylabel=labelsolte[2],
    fontfamily="Times",linecolor=c1,framestyle=:box,legend=:topright,
    size=plot_size,margin=plot_margin)
plot!(ers,XFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 0.3),label="",
    linecolor=c1)
plot!(ers,XFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 0.3),
    label=L"M(t)",linecolor=c1)
plot!(ers,QFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 0.3),label="",
    linecolor=c2)
plot!(ers,QFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 0.3),label="",
    linecolor=c2)
plot!(ers,QFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 0.3),
    label=L"Q(t)",linecolor=c2)

savefig(dir*"TimeSoln_HW_3.pdf")

#-----------------------------------------------------------------------
#   FIGURE 4
#-----------------------------------------------------------------------

σM  = 0.3*ONE

p   = RSHWgen(tpfl(1.0),tpfl(1.0),σM,tpfl,0)

s0  = p[1]
Ms0 = σM*(1e8)
Ms  = tpfl(1474.07)*Ms0

u0  = [ tpfl(1.0) , tpfl(0.1) ]

f   = (u,p,t)->HWslv.F(u,p)

tspan1 = tpfl.((0.0,1.5*Δτest(u0[1],p)))  ;

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

c1=:steelblue3
c2=:red3
c3=:darkorange3

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

plot(X1s,Y1s,tpfl(0),tpfl(1e-2),xlims=(0,1.02),ylims=ylms,
    label=L"Y(\mu)",
    xlabel=labelphase[1],ylabel=labelphase[2],fontfamily="Times",
    linecolor=c1,framestyle=:box,legend=:topleft,size=plot_size,
    margin=plot_margin)
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)

savefig(dir*"PhaseSpace_HW_4.pdf")

    TF          = tf1
    TFyr        = Float64(TF*Ms/s0)*(1/(3e8))*(1/3.154e7)

    XFFn    = t->Xr1(TF*t)
    YFFn    = t->Yr1(TF*t)

    XFFns   = s->XFFn(ers(s))
    YFFns   = s->YFFn(ers(s))

    QFFn    = t->sqrt(YFFn(t))*XFFn(t)
    QFFns   = s->sqrt(YFFns(s))*XFFns(s)

    labelsolte = ( snst(TFyr) , snsMu(Ms0,2) )

plot(ers,XFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 1),label="",
    xlabel=labelsolte[1],ylabel=labelsolte[2],
    fontfamily="Times",linecolor=c1,framestyle=:box,legend=:topright,
    size=plot_size,margin=plot_margin)
plot!(ers,XFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c1)
plot!(ers,XFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 1),
    label=L"M(t)",linecolor=c1)
plot!(ers,QFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c2)
plot!(ers,QFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c2)
plot!(ers,QFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 1),label=L"Q(t)",
    linecolor=c2)

savefig(dir*"TimeSoln_HW_4.pdf")

#-----------------------------------------------------------------------
#   FIGURE 5
#-----------------------------------------------------------------------

σM  = 0.9*ONE

p   = RSHWgen(tpfl(1.0),tpfl(1.0),σM,tpfl,0)

s0  = p[1]
Ms0 = σM*(1e8)
Ms  = tpfl(1474.07)*Ms0

u0  = [ tpfl(1.0) , tpfl(0.1) ]

f   = (u,p,t)->HWslv.F(u,p)

tspan1 = tpfl.((0.0,1.5*Δτest(u0[1],p)))

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

c1=:steelblue3
c2=:red3
c3=:darkorange3

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

plot(X1s,Y1s,tpfl(0),tpfl(1e-2),xlims=(0,1.02),ylims=ylms,
    label=L"Y(\mu)",xlabel=labelphase[1],ylabel=labelphase[2],
    fontfamily="Times",linecolor=c1,framestyle=:box,
    legend=(0.25,0.25),size=plot_size,margin=plot_margin)
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)

savefig(dir*"PhaseSpace_HW_5.pdf")

    TF          = tf1
    TFyr        = Float64(TF*Ms/s0)*(1/(3e8))*(1/3.154e7)

    XFFn    = t->Xr1(TF*t)
    YFFn    = t->Yr1(TF*t)

    XFFns   = s->XFFn(ers(s))
    YFFns   = s->YFFn(ers(s))

    QFFn    = t->sqrt(YFFn(t))*XFFn(t)
    QFFns   = s->sqrt(YFFns(s))*XFFns(s)

    labelsolte = ( snst(TFyr) , snsMu(Ms0,2) )

plot(ers,XFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 1),label="",
    xlabel=labelsolte[1],ylabel=labelsolte[2],
    fontfamily="Times",linecolor=c1,framestyle=:box,legend=:topright,
    size=plot_size,margin=plot_margin)
plot!(ers,XFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c1)
plot!(ers,XFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 1),
    label=L"M(t)",linecolor=c1)
plot!(ers,QFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c2)
plot!(ers,QFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c2)
plot!(ers,QFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 1),label=L"Q(t)",
    linecolor=c2)

savefig(dir*"TimeSoln_HW_5.pdf")

#-----------------------------------------------------------------------
#   FIGURE 6
#-----------------------------------------------------------------------

σM  = 1.5*ONE

p   = RSHWgen(tpfl(1.0),tpfl(1.0),σM,tpfl,0)

s0  = p[1]
Ms0 = σM*(1e8)
Ms  = tpfl(1474.07)*Ms0

u0  = [ tpfl(1.0) , tpfl(0.1) ]

f   = (u,p,t)->HWslv.F(u,p)

tspan1 = tpfl.((0.0,1.5*Δτest(u0[1],p)))

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

c1=:steelblue3
c2=:red3
c3=:darkorange3

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

plot(X1s,Y1s,tpfl(0),tpfl(1e-2),xlims=(0,1.02),ylims=ylms,
    label=L"Y(\mu)",xlabel=labelphase[1],ylabel=labelphase[2],
    fontfamily="Times",linecolor=c1,framestyle=:box,legend=(0.25,0.25),
    size=plot_size,margin=plot_margin)
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)

savefig(dir*"PhaseSpace_HW_6.pdf")

    TF          = tf1
    TFyr        = Float64(TF*Ms/s0)*(1/(3e8))*(1/3.154e7)

    XFFn    = t->Xr1(TF*t)
    YFFn    = t->Yr1(TF*t)

    XFFns   = s->XFFn(ers(s))
    YFFns   = s->YFFn(ers(s))

    QFFn    = t->sqrt(YFFn(t))*XFFn(t)
    QFFns   = s->sqrt(YFFns(s))*XFFns(s)

    labelsolte = ( snst(TFyr) , snsMu(Ms0,2) )

plot(ers,XFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 1),label="",
    xlabel=labelsolte[1],ylabel=labelsolte[2],
    fontfamily="Times",linecolor=c1,framestyle=:box,legend=:topright,
    size=plot_size,margin=plot_margin)
plot!(ers,XFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c1)
plot!(ers,XFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 1),
    label=L"M(t)",linecolor=c1)
plot!(ers,QFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c2)
plot!(ers,QFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c2)
plot!(ers,QFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 1),label=L"Q(t)",
    linecolor=c2)

savefig(dir*"TimeSoln_HW_6.pdf")

#-----------------------------------------------------------------------
#   FIGURE 7
#-----------------------------------------------------------------------

σM  = 1.68*ONE

p   = RSHWgen(tpfl(1.0),tpfl(1.0),σM,tpfl,0)

s0  = p[1]
Ms0 = σM*(1e8)
Ms  = tpfl(1474.07)*Ms0

u0  = [ tpfl(1.0) , tpfl(0.1) ]

f   = (u,p,t)->HWslv.F(u,p)

tspan1 = tpfl.((0.0,1.5*Δτest(u0[1],p)))

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

c1=:steelblue3
c2=:red3
c3=:darkorange3

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

plot(X1s,Y1s,tpfl(0),tpfl(6e-3),xlims=(0,1.02),ylims=ylms,
    label=L"Y(\mu)",xlabel=labelphase[1],ylabel=labelphase[2],
    fontfamily="Times",linecolor=c1,framestyle=:box,
    legend=(0.25,0.25),size=plot_size,margin=plot_margin)
plot!(X1s,Y1s,tpfl(6e-3),tpfl(6.3e-3),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(6.3e-3),tpfl(6.4e-3),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(6.4e-3),tpfl(2e-2),xlims=(0,1.02),ylims=ylms,
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)
    
savefig(dir*"PhaseSpace_HW_7.pdf")

    TF          = tf1
    TFyr        = Float64(TF*Ms/s0)*(1/(3e8))*(1/3.154e7)

    XFFn    = t->Xr1(TF*t)
    YFFn    = t->Yr1(TF*t)

    XFFns   = s->XFFn(ers(s))
    YFFns   = s->YFFn(ers(s))

    QFFn    = t->sqrt(YFFn(t))*XFFn(t)
    QFFns   = s->sqrt(YFFns(s))*XFFns(s)

    labelsolte = ( snst(TFyr) , snsMu(Ms0,2) )

    plot(ers,XFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 1),label="",
            xlabel=labelsolte[1],
            ylabel=labelsolte[2],
            fontfamily="Times",linecolor=c1,framestyle=:box,
            legend=:topright,size=plot_size,margin=plot_margin)
    plot!(ers,XFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 1),label="",
          linecolor=c1)
    plot!(ers,XFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 1),
          label=L"M(t)",linecolor=c1)
    plot!(ers,QFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 1),label="",
          linecolor=c2)
    plot!(ers,QFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 1),label="",
          linecolor=c2)
    plot!(ers,QFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 1),
          label=L"Q(t)",linecolor=c2)

    savefig(dir*"TimeSoln_HW_7.pdf")

#-----------------------------------------------------------------------
#   FIGURE 8
#-----------------------------------------------------------------------

σM  = 1.0*ONE

p   = RSHWgen(tpfl(1.0),tpfl(1.0),σM,tpfl,0)

s0  = p[1]
Ms0 = σM*(1e8)
Ms  = tpfl(1474.07)*Ms0

u0  = [ tpfl(1.0) , tpfl(0.7)^2 ]

f   = (u,p,t)->HWslv.F(u,p)

tspan1 = tpfl.((0.0,1.5*Δτest(u0[1],p)))

sol = solve(   ODEProblem(f,u0,tspan1,p),integrator,
                reltol=rt,abstol=at,callback=cb)

tf1  = sol.t[end]

Xr1  = t->(sol(t)[1])
Yr1  = t->(sol(t)[2])

X1   = t->Xr1(t*tf1)
Y1   = t->Yr1(t*tf1)

X1s  = s->X1(ers(s))
Y1s  = s->Y1(ers(s))

ylms=(0,1.05)

c1=:steelblue3
c2=:red3
c3=:darkorange3

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

plot(s->0.7,s->s,tpfl(0),tpfl(0.995),xlims=(0,1.02),ylims=ylms,
    label=L"\mu_h",linestyle=:dash,linecolor=:orange,
    xlabel=L"\mu=M/(1.0\times 10^{8}M_\odot)",ylabel=L"Y=(Q/M)^2",
    fontfamily="Times",framestyle=:box,legend=(0.30,0.40),
    size=plot_size,margin=plot_margin)
plot!(s->p[2],s->s,tpfl(0),tpfl(0.997),xlims=(0,1.02),ylims=ylms,
    label=L"z_0",linestyle=:dot,linecolor=:grey50)
plot!(μ->HWslv.YH(μ,(0.7)^2),tpfl(0.7),tpfl(1),linestyle=:dashdot,
    label=L"Y_H(\mu)",linecolor=:grey20)
plot!(μ->tpfl(1),p[2],tpfl(0.7),xlims=(0, 1),ylims=ylms,label="",
    linestyle=:dashdot,linecolor=:grey20)
plot!(Y->HWslv.muattr(Y,p),Y->Y,tpfl(0),tpfl(0.98),xlims=(0, 1),
    ylims=ylms,label="",linestyle=:dashdot,linecolor=:grey20)
plot!(Y->HWslv.muattr(Y,p),Y->Y,tpfl(0.98),tpfl(1),xlims=(0, 1),
    ylims=ylms,label="",linestyle=:dashdot,linecolor=:grey20)
plot!(X1s,Y1s,tpfl(0),tpfl(6e-3),xlims=(0,1.02),ylims=ylms,
    label=L"\mathrm{Num.}",linecolor=c1)
plot!(X1s,Y1s,tpfl(6e-3),tpfl(6.3e-3),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(6.3e-3),tpfl(6.4e-3),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(6.4e-3),tpfl(2e-2),xlims=(0,1.02),ylims=ylms,
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)

savefig(dir*"PhaseSpace_HW_8.pdf")

    TF          = tf1
    TFyr        = Float64(TF*Ms/s0)*(1/(3e8))*(1/3.154e7)

    XFFn    = t->Xr1(TF*t)
    YFFn    = t->Yr1(TF*t)

    XFFns   = s->XFFn(ers(s))
    YFFns   = s->YFFn(ers(s))

    QFFn    = t->sqrt(YFFn(t))*XFFn(t)
    QFFns   = s->sqrt(YFFns(s))*XFFns(s)

    labelsolte = ( snst(TFyr) , snsMu(Ms0,2) )

plot(ers,XFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 1),label="",
    xlabel=labelsolte[1],ylabel=labelsolte[2],
    fontfamily="Times",linecolor=c1,framestyle=:box,legend=:topright,
    size=plot_size,margin=plot_margin)
plot!(ers,XFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c1)
plot!(ers,XFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 1),
    label=L"M(t)",linecolor=c1)
plot!(ers,QFFns,0,0.01,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c2)
plot!(ers,QFFns,0.01,0.98,xlims=(0, 1.02),ylims=(0, 1),label="",
    linecolor=c2)
plot!(ers,QFFns,0.98,0.9999,xlims=(0, 1.02),ylims=(0, 1),label=L"Q(t)",
    linecolor=c2)

savefig(dir*"TimeSoln_HW_8.pdf")

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

function muanalytical(s, p)
    Y = Y1s(s)
    return HWslv.muattr(Y, p)
end

residual8 = s->abs(X1s(s) - muanalytical(s, p))

residual8rl = s->abs(X1s(s) - muanalytical(s, p)) / 
                abs(muanalytical(s, p))

plot(s->Y1s(s),residual8,tpfl(1e-2), tpfl(0.9938),
    xlabel=L"Y=(Q/M)^2", ylabel=L"\Delta \mu",
    label="",grid=true,legend=false,
    fontfamily="Times",linecolor=c1,framestyle=:box,
    size=plot_size,margin=plot_margin,xlims=(0, 1), ylims=(0.0,0.004))
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.9938), 
        tpfl(0.99395),linecolor=c1)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.99395), 
        tpfl(0.99409),linecolor=c1)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.99409),
        tpfl(0.99415),linecolor=c1)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.99415),
        tpfl(1)-tpfl(3.5e-3),linecolor=c1)

savefig(dir*"PhaseSpace_HW_8_residual.pdf")

plot(s->Y1s(s),residual8,tpfl(1e-2), tpfl(0.9938),
    xlabel=L"Y=(Q/M)^2", ylabel=L"\Delta \mu",
    label="",grid=true,legend=false,
    fontfamily="Times",linecolor=c1,framestyle=:box,
    size=plot_size,margin=plot_margin,xlims=(1e-4, 1), 
    ylims=(1e-6,0.1),yscale=:log10)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.9938), 
        tpfl(0.99395),linecolor=c1)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.99395), 
        tpfl(0.99409),linecolor=c1)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.99409), 
        tpfl(0.99415),linecolor=c1)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.99415),
        tpfl(1)-tpfl(3.5e-3),linecolor=c1)

savefig(dir*"PhaseSpace_HW_8_residual_log.pdf")

plot(s->Y1s(s),residual8rl,tpfl(1e-2), tpfl(0.9938),
    xlabel=L"Y=(Q/M)^2", ylabel=L"\Delta \mu",
    label="",grid=true,legend=false,
    fontfamily="Times",linecolor=c1,framestyle=:box,
    size=plot_size,margin=plot_margin,xlims=(1e-4, 1), 
    ylims=(1e-6,0.1),yscale=:log10)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.9938), 
        tpfl(0.99395),linecolor=c1)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.99395), 
        tpfl(0.99409),linecolor=c1)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.99409), 
        tpfl(0.99415),linecolor=c1)
plot!(s->Y1s(s), s->abs(X1s(s)-muanalytical(s,p)), tpfl(0.99415),
        tpfl(1)-tpfl(3.5e-3),linecolor=c1)

savefig(dir*"PhaseSpace_HW_8_residual_log_rel.pdf")

#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
#
#   CONFIGURATION SPACE PLOTS
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

setprecision(ArbFloat, 400)
tpfl    = ArbFloat{400}

plot_margin = 3Plots.mm
plot_size = (340, 240)

ONE = one(tpfl)

σM  = tpfl(1.0*ONE)

p   = RSHWgen(tpfl(1.0),tpfl(1.0),σM,tpfl)

s0  = p[1]
Ms0 = σM*(1e8)
Ms  = tpfl(1474.07)*Ms0

f   = (u,p,t)->HWslv.F(u,p)

c1=:royalblue3
c2=:red3
c3=:darkorange3

#-----------------------------------------------------------------------

u0  = [ tpfl(0.15) , tpfl(0.998) ]
tspan1 = tpfl.((0.0,ONE/2))

sol = solve(   ODEProblem(f,u0,tspan1,p),integrator,
                reltol=rt,abstol=rt,callback=cb)

tf1  = sol.t[end]

Xr1  = t->(sol(t)[1])
Yr1  = t->(sol(t)[2])

X1   = t->Xr1(t*tf1)
Y1   = t->Yr1(t*tf1)

X1s  = s->X1(ers(s))
Y1s  = s->Y1(ers(s))

ylms=(0,1.04)

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

plot(X1s,Y1s,tpfl(0),tpfl(1e-2),xlims=(0,1.02),ylims=ylms,
    label="",
    xlabel=labelphase[1],
    ylabel=labelphase[2],
    fontfamily="Times",
    linecolor=c1,
    framestyle=:box,
    size=plot_size,
    margin=plot_margin)
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)

#-----------------------------------------------------------------------

u0  = [ tpfl(0.27) , tpfl(0.998) ]
tspan1 = tpfl.((0.0,ONE/2))

sol = solve(   ODEProblem(f,u0,tspan1,p),integrator,
                reltol=rt,abstol=rt,callback=cb)

tf1  = sol.t[end]

Xr1  = t->(sol(t)[1])
Yr1  = t->(sol(t)[2])

X1   = t->Xr1(t*tf1)
Y1   = t->Yr1(t*tf1)

X1s  = s->X1(ers(s))
Y1s  = s->Y1(ers(s))

ylms=(0,1.04)

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

c1=:purple3

plot!(X1s,Y1s,tpfl(0),tpfl(1e-2),xlims=(0,1.02),ylims=ylms,
    label="",
    xlabel=labelphase[1],
    ylabel=labelphase[2],
    fontfamily="Times",
    linecolor=c1,
    framestyle=:box,
    size=plot_size,
    margin=plot_margin)
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)

#-----------------------------------------------------------------------

setprecision(ArbFloat, 200)
tpfl    = ArbFloat{200}

ONE = one(tpfl)

σM  = tpfl(1.0*ONE)

p   = RSHWgen(tpfl(1.0),tpfl(1.0),σM,tpfl)

s0  = p[1]
Ms0 = σM*(1e8)
Ms  = tpfl(1474.07)*Ms0

f   = (u,p,t)->HWslv.F(u,p)

#-----------------------------------------------------------------------

u0  = [ tpfl(1.0) , tpfl(0.04) ]
tspan1 = tpfl.((0.0,1.5*Δτest(u0[1],p)))

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

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

c1=:darkred

plot!(X1s,Y1s,tpfl(0),tpfl(1e-2),xlims=(0,1.02),ylims=ylms,
    label="",
    xlabel=labelphase[1],
    ylabel=labelphase[2],
    fontfamily="Times",
    linecolor=c1,
    framestyle=:box,
    size=plot_size,
    margin=plot_margin)
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)

#-----------------------------------------------------------------------

u0  = [ tpfl(1.0) , tpfl(0.2) ]
tspan1 = tpfl.((0.0,1.5*Δτest(u0[1],p)))

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

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

c1=:darkorange3

plot!(X1s,Y1s,tpfl(0),tpfl(1e-2),xlims=(0,1.02),ylims=ylms,
    label="",
    xlabel=labelphase[1],
    ylabel=labelphase[2],
    fontfamily="Times",
    linecolor=c1,
    framestyle=:box,
    size=plot_size,
    margin=plot_margin)
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)

#-----------------------------------------------------------------------

u0  = [ tpfl(1.0) , tpfl(0.4) ]
tspan1 = tpfl.((0.0,1.5*Δτest(u0[1],p)))

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

labelphase = ( snsM(Ms0,2) , L"(Q/M)^2" )

c1=:darkgreen

plot!(X1s,Y1s,tpfl(0),tpfl(1e-2),xlims=(0,1.02),ylims=ylms,
    label="",
    xlabel=labelphase[1],
    ylabel=labelphase[2],
    fontfamily="Times",
    linecolor=c1,
    framestyle=:box,
    size=plot_size,
    margin=plot_margin)
plot!(X1s,Y1s,tpfl(1e-2),tpfl(1.05e-2),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1.05e-2),tpfl(1.1e-2),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1.1e-2),tpfl(1.6e-2),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1.6e-2),tpfl(2e-2),xlims=(0,1.02),ylims=ylms,
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
    label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(0.87),tpfl(1)-tpfl(8e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(8e-2),tpfl(1)-tpfl(6e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(6e-2),tpfl(1)-tpfl(5e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(5e-2),tpfl(1)-tpfl(4e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(4e-2),tpfl(1)-tpfl(3.1e-2),xlims=(0,1.02),
    ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.1e-2),tpfl(1)-tpfl(1.8e-2),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(1.8e-2),tpfl(1)-tpfl(7.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(7.5e-3),tpfl(1)-tpfl(3.5e-3),
    xlims=(0,1.02),ylims=ylms,label="",linecolor=c1)
plot!(X1s,Y1s,tpfl(1)-tpfl(3.5e-3),tpfl(1),xlims=(0,1.02),ylims=ylms,
    label="",linecolor=c1)

#-----------------------------------------------------------------------

savefig(dir*"PhaseSpace_HW_Mult.pdf")

#-----------------------------------------------------------------------