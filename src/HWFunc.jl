#-----------------------------------------------------------------------
#
#   GENERAL FUNCTIONS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    integral_ode( f::Function, tspan, p; abstol=1e-10, reltol=1e-8, 
                                         solver=Tsit5() )

This function computes an integral using a numerical ode solver (instead
of using quadrature methods).
"""
function integral_ode( f::Function, tspan, p; abstol=1e-12, reltol=1e-12, 
                                              solver=Vern9() )
    function intgrnd!(du, u, p, x)
        du[1] = f(x) 
    end
    sol = solve( ODEProblem(intgrnd!,[tpfl(0.0)],tspan,p) , solver; 
                 abstol=abstol, reltol=reltol , save_everystep=false )
    return sol.u[end][1]
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    hvs(x)

This is the Heaviside function. Returns '1' for positive arguments, '0' 
for negative arguments and '1/2' if the argument is zero.
"""
function hvs(x)
    tpfl = typeof(x)
    (ONE,ZER,HALF)  = (one(tpfl),zero(tpfl),tpfl(1/2))
    θ = (x < ZER) ? ZER : ONE
    return (x == ZER) ? HALF : θ
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    ers(x)

This is the Heaviside function. Returns '1' for positive arguments, '0' 
for negative arguments and '1/2' if the argument is zero.
"""            #   ers=x->(1 + erf((3*(2*x-1))/(4*sqrt(-((x-1)*x)))))/2;
function ers(x)
    tpfl = typeof(x)
    (ONE,TWO,THR,FOU) = (one(tpfl),tpfl(2),tpfl(3),tpfl(4))
    return (ONE + erf((THR*(TWO*x-ONE))/(FOU*sqrt(-((x-ONE)*x)))))/TWO;
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   PARAMETER FUNCTIONS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    pars(tpfl=Float64,αn=-1)
        
This function sets the values of the parameters `s0`, `z0` and `b0`,
as well as the constants used throughout. Unit conventions may be found
in appendix C of (arXiv:2503.20696).
"""
function pars(tpfl=Float64,αn=-1)
    # Numbers
    (ONE,TWO,THR,FOU,FIV) = (one(tpfl),tpfl(2),tpfl(3),tpfl(4),tpfl(5))
    (SIX,SEV,EIG,NIN,PI)  = (tpfl(6),tpfl(7),tpfl(8),tpfl(9),tpfl(π))

    α1    = tpfl(0.26792)                 # No massless neutrinos
    α2    = tpfl(2.0228)                  # Three massless neutrinos

    ħ     = tpfl(2.60624e-70)             # Planck constant (m^2)
    qe    = tpfl(1.38452e-36)             # Electron charge (m)
    me    = tpfl(6.75107e-58)             # Electron mass (m)

    Msol  = tpfl(1477)                    # Solar mass (meters)
    Ms    = tpfl(1e8)*Msol                # Mass scale
    tuniv = tpfl(1.3043241268684672e26)   # Age of universe (meters)

    c     = tpfl(299792458)               # Speed of light
    yr2s  = tpfl(3.15576e7)               # Year to second conversion
    yr2m  = tpfl(9.460536207068014e15)    # Year to meter conversion

    if αn < 0
        α = α1
    elseif αn == 0 
        α = α2
    else
        α =  αn
    end
    
    s0    = (α*ħ)/( tpfl(1920)*PI*Ms^2 ) 
    z0    = ((qe*ħ)/(Ms*PI*me^2)) * log(
                    tpfl(960)*(qe^4)*(Ms^2)/(((me^2)*(PI^2)*α*(ħ^2)))
                    )
    b0    = (me^2)*Ms*PI/(qe*ħ)

    return ( (s0,z0,b0) , (qe,me,ħ) , (Ms,Msol,tuniv) , (α1,α2) ,
             (ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI) , (c,yr2s,yr2m) )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   RESCALING FUNCTIONS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    ξF(σm,σet,tpfl=Float64)
    
This function computes the charge to mass scaling parameter 'ξ'.
"""
function ξF(σm,σet,tpfl=Float64,αn=-1)
    (s0,z0,b0)    = pars(tpfl,αn)[1]
    (ZER,ONE,TWO) = (zero(tpfl),one(tpfl),tpfl(2))
    if σm>ZER && σet>ZER
        ϑ=(σet/σm)^3
        den = b0*z0+TWO*log(ϑ)
        if den<=ZER
            return -ONE
        else
            return (b0*z0*ϑ/den)^(1/3)
        end
    else
        return -ONE
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    σMF(σm,σet,tpfl=Float64)
    
This function computes the adjustment 'σM' to the mass scale 'Ms'.
"""
function σMF(σm,σet,tpfl=Float64,αn=-1)
    (ZER,ONE) = (zero(tpfl),one(tpfl))
    if σm>ZER && σet>ZER
        ϑ=(σet/σm)^3
        ξ = ξF(σm,σet,tpfl,αn)
        if ξ<=ZER
            return -ONE
        else
            return ϑ/(σm*ξ^2)
        end
    else
        return -ONE
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    sigtransF(σm,σet,tpfl=Float64)
    
This function computes the true charge rescaling parameter `σe` from the
charge pseudoscaling parameter `σet` and the dimensionless threshold 
mass parameter `σm`.
"""
function sigtransF(σm,σet,tpfl=Float64)
    (ZER,ONE) = (zero(tpfl),one(tpfl))
    ξ   = ξF(σm,σet,tpfl)
    if ξ<=ZER
        return (-ONE,-ONE)
    else
        return (σm,σm*ξ)
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    RSHW(σm,σet,tpfl=Float64,αn=-1)
    
This function computes the rescaling of the parameters `s0`, `b0` and
`σM` under a rescaling that preserves `z0`.
"""
function RSHW(σm,σet,tpfl=Float64,αn=-1)
    (s0,z0,b0)    = pars(tpfl,αn)[1]
    (ZER,ONE,TWO) = (zero(tpfl),one(tpfl),tpfl(2))
    if σm>ZER && σet>ZER
        ϑ = (σet/σm)^3
        ξ = ξF(σm,σet,tpfl,αn)
        if ξ<=ZER
            return (s0,z0,b0,ONE)
        else
            σM=σMF(σm,σet,tpfl,αn)
            return ( s0/σM^2 , z0 , b0+TWO*log(ϑ)/z0 , σM )
        end
    else
        return (s0,z0,b0,ONE)
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    RSHWgen(σm,σe,σM,tpfl=Float64)
    
This function computes the rescaling of the parameters `s0`, `z0`, and 
`b0` under a general rescaling given by the rescaling parameters `σm`,
`σe` and `σM`.
"""
function RSHWgen(σm,σe,σM,tpfl=Float64,αn=-1)
    (s0,z0,b0)    = pars(tpfl,αn)[1]
    (ZER,ONE,TWO) = (zero(tpfl),one(tpfl),tpfl(2))
    if σm>ZER && σe>ZER && σM>ZER
        Ξ = ( z0 + (TWO/b0)*log(((σe^2)*σM)/σm) )
        if Ξ<=ZER
            return (s0,z0,b0)
        else
            S0  = s0/(σM^2) 
            Z0  = (σe/((σm^2)*σM)) * Ξ
            B0  = (σm^2)*σM*b0 / σe
            return ( S0 , Z0 , B0 )
        end
    else
        return (s0,z0,b0)
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   TIME ESTIMATE FUNCTIONS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Integrand for Δτatt
#-----------------------------------------------------------------------
"""
    dτhwk(μ,μ1,z0)

This function supplies the integrand for 'Δτ=∫dτ' on the mass 
dissipation curve, for an initial charge to mass ratio 'μ1'.
"""
function dτhwk(μ,μ1,z0)
    tpfl = typeof(μ)
    (ONE,THR,FOU,EIG,NIN) = (tpfl(1),tpfl(3),tpfl(4),tpfl(8),tpfl(9))
    (FRT,SXT) = (tpfl(14),tpfl(16))
    Y1      = μ1^2
    (A1,A2) = (abs(μ^2-Y1),abs(NIN*μ^2-EIG*Y1))
    return  ( (Y1^2 + FOU*sqrt(A1)*(TWO*A1+Y1)*μ + EIG*A1*μ^2)^2 * 
              (TWO*A1 + μ*(sqrt(A2) + μ)) 
            )/(A1^2*(sqrt(A2) + THR*μ)^4)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Δτhwk Integral 
#-----------------------------------------------------------------------
"""
    Δτhwk(μ1,z0)

The integral 'Δτ=∫dτ' from one to the attractor curve for a given value 
of 'z0'. This estimates the amount of rescaled time 'τ' that the system 
spends on the mass dissipation curve, for an initial charge to mass 
ratio 'μ1'.
"""
function Δτhwk(μ1,z0)
    tpfl    = typeof(μ1)
    (ONE,TWO,FOU) = (tpfl(1),tpfl(2),tpfl(4))
    HALF = ONE/TWO
    if μ1<z0 
        x = sqrt( HALF+sqrt(z0*μ1) - sqrt(ONE-FOU*μ1^2+FOU*sqrt(z0*μ1))/
                  TWO)
        return integral_ode(y -> dτhwk(y,μ1,z0),x,ONE)
    elseif μ1>z0
        return integral_ode(y -> dτhwk(y,μ1,z0),μ1,ONE)
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Integrand for Δτatt
#-----------------------------------------------------------------------
"""
    dτatt(x, z0)

This function supplies the integrand for 'Δτ=∫dτ' on the attractor 
curve.
"""
function dτatt(x, z0)
    tpfl = typeof(x)
    (ONE,THR,FOU,EIG,NIN) = (tpfl(1),tpfl(3),tpfl(4),tpfl(8),tpfl(9))
    (FRT,SXT) = (tpfl(14),tpfl(16))
    if x == zero(tpfl)
        return zero(tpfl)
    end
    P = THR*(ONE+x) + sqrt(-FRT*x+NIN*(ONE+x^2))
    return -( ( sqrt(x) / (FOU*(ONE+THR*x)*z0) ) * 
              ( (P^4*(ONE-x)^4) / (FOU*x*(EIG*x - P*(ONE+x))*z0^2) - 
                SXT*x^2
              )
            )^(-1)
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Δτatt Integral 
#-----------------------------------------------------------------------
"""
    Δτatt(μ,z0)

The integral 'Δτ=∫dτ' from zero to 'μ' for a given value of 'z0'. This
estimates the amount of rescaled time 'τ' that the system spends on the
attractor curve.
"""
function Δτatt(μ,z0)
    tpfl    = typeof(μ)
    ETO     = tpfl(81)
    (ONE,TWO,THR,SIX,NIN) = (tpfl(1),tpfl(2),tpfl(3),tpfl(6),tpfl(9))
    (n1,n2) = (ONE/THR,TWO/THR)
    x = (-z0^n2 + (z0^2 + SIX*μ*(NIN*μ + sqrt(THR*z0^2+ETO*μ^2)))^n1)^2/
        (THR*(z0^n2)*(z0^2 + SIX*μ*(NIN*μ +sqrt(THR*z0^2+ETO*μ^2)))^n1)
    return integral_ode(y -> dτatt(y, z0),x,tpfl(0.0))
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Δτatt Integral for small μ
#-----------------------------------------------------------------------
"""
    Δτatts(μ)

The integral 'Δτ=∫dτ' from zero to 'μ', valid for 'μ≪1'.
"""
function Δτatts(μ)
    return (tpfl(32)/tpfl(81)) * μ^3
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Evaporation rescaled time estimate
#-----------------------------------------------------------------------
"""
    Δτfull(μ,p)

This function estimates the evaporation timescales for an initial 
charge to mass ratio 'μ' and a mass of 'M=Ms'.
"""
function Δτfull(μ,p)
    tpfl = typeof(μ)
    (s0,z0,b0) = p
    if μ >= z0
        return ((exp(b0*(μ-z0))-one(tpfl))/b0+Δτatt(z0,z0)+Δτhwk(μ,z0))
    elseif μ >= tpfl(1e-4)
        return (Δτatt(μ,z0)+Δτhwk(μ,z0))
    elseif μ >= zero(tpfl)
        return (Δτatts(μ)+Δτhwk(μ,z0))
    else
        return zero(tpfl)
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Evaporation rescaled time estimate
#-----------------------------------------------------------------------
"""
    Δτane(μ,p)

This function estimates the evaporation timescales for an initial
rescaled mass 'μ', assuming the system initially is either near-extremal
or close to the approximate attractor curve.
"""
function Δτane(μ,p)
    tpfl = typeof(μ)
    (s0,z0,b0) = p
    if μ >= z0
        return ((exp(b0*(μ-z0))-one(tpfl))/b0+Δτatt(z0,z0))
    elseif μ >= tpfl(1e-4)
        return (Ms/(s0*yr2m)) * (Δτatt(μ,z0))
    elseif μ >= zero(tpfl)
        return (Ms/(s0*yr2m)) * (Δτatts(μ))
    else
        return zero(tpfl)
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Evaporation rescaled time estimate
#-----------------------------------------------------------------------
"""
    Δτest(μ,p)

This function provides a simple estimate for the evaporation timescales 
for an initial charge to mass ratio 'μ' and parameters 'p'.
"""
function Δτest(μ,p)
    tpfl = typeof(μ)
    (s0,z0,b0) = p
    if μ >= z0
        return ((exp(b0*(μ-z0))-one(tpfl))/b0+Δτatt(z0,z0))
    else
        return Δτatt(z0,z0)
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Evaporation time estimate
#-----------------------------------------------------------------------
"""
    Δt(μ,σm,σet,tpfl=Float64)

This function estimates the evaporation timescales for an initial 
rescaled mass 'μ', assuming the system initially is either near-extremal 
or close to the approximate attractor curve.
"""
function Δt(μ,σm,σet,tpfl=Float64,αn=-1)
    p    = pars(tpfl,αn) 
    (Ms0,Msol,tuniv)  = p[3]
    (c,yr2s,yr2m)     = p[6]
    (s0,z0,b0,σM) = RSHW(σm,σet,tpfl,αn)
    Ms = σM*Ms0
    if μ >= z0
        Δτ0 = Δτatt(z0,z0)
        return (Ms/(s0*yr2m)) * ((exp(b0*(μ-z0))-one(tpfl))/b0 + Δτ0)
    elseif μ >= tpfl(1e-4)
        return (Ms/(s0*yr2m)) * Δτatt(μ,z0)
    elseif μ >= zero(tpfl)
        return (Ms/(s0*yr2m)) * Δτatts(μ)
    else
        return zero(tpfl)
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Simple evaporation time estimate
#-----------------------------------------------------------------------
"""
    Δts(μ,σm,σet,tpfl=Float64)

This function estimates the evaporation timescales for an initial 
rescaled mass 'μ', assuming the system initially is either near-extremal 
with 'μ-z0≫1/|b0|'.
"""
function Δts(μ,σm,σet,tpfl=Float64,αn=-1)
    p    = pars(tpfl,αn) 
    (Ms0,Msol,tuniv)  = p[3]
    (c,yr2s,yr2m)     = p[6]
    (s0,z0,b0,σM) = RSHW(σm,σet,tpfl,αn)
    Ms = σM*Ms0
    if μ >= z0 && abs(b0*(μ-z0))>one(tpfl)
        Δτ0 = Δτatt(z0,z0)
        return (Ms/(s0*yr2m)) * (exp(b0*(μ-z0)))/b0
    else
        return zero(tpfl)
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#   Universe age to rescaled time
#-----------------------------------------------------------------------
"""
    τuniv(σm,σet,tpfl=Float64)

This function computes the age of the universe in rescaled time units.
"""
function τuniv(σm,σet,tpfl=Float64,αn=-1)
    p    = pars(tpfl,αn) 
    (Ms0,Msol,tuniv)  = p[3]
    (c,yr2s,yr2m)     = p[6]
    if σm>ZER && σet>ZER
        ξ = ξF(σm,σet,tpfl)
        if ξ<=ZER
            return ZER
        else
            (s0,z0,b0,σM) = RSHW(σm,σet,tpfl,αn)
            return s0*tuniv/(σM*Ms0)
        end
    else
        return ZER
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
#
#   MASS RESCALING FUNCTIONS
#
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    μFromMass(M,σm,σet,tpfl=Float64)
    
This function computes the rescaled mass 'μ' from the mass 'M' for a
given set of charge pseudoscaling parameters 'σm' and 'σet'.
"""
function μFromMass(M,σm,σet,tpfl=Float64,αn=-1)
    p = pars(tpfl,αn)
    (Ms,Msol,tuniv) = p[3]
    (ZER,ONE,TWO) = (zero(tpfl),one(tpfl),tpfl(2))
    if σm>ZER && σet>ZER
        ξ = ξF(σm,σet,tpfl,αn)
        if ξ<=ZER
            return M/Ms
        else
            return M/(σMF(σm,σet,tpfl,αn)*Ms)
        end
    else
        return M/Ms
    end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   μ from τ
#-----------------------------------------------------------------------
"""
    μfromτ(τ,b0,z0=tpfl(0.5313221978739817),
            Δτattz0=tpfl(0.30505917474345556),tpfl=Float64)

This function computes the rescaled mass 'μ' from the rescaled time 'τ'.
"""
function μfromτ(τ,b0,z0=tpfl(0.5313221978739817),
                    Δτattz0=tpfl(0.30505917474345556),tpfl=Float64)
    ONE     = one(tpfl)
    if τ <= Δτattz0
        return zero(tpfl)
    elseif τ > Δτattz0
        return z0+log(ONE+b0*(τ-Δτattz0))/b0
    else
        return zero(tpfl)
    end
end #-------------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    Mz0( σm , σet , tpfl=Float64 , αn=0 )
    
This function computes the mass associated with a given value of the
dimensionless threshold mass parameter `z0`.
"""
function Mz0( σm , σet , tpfl=Float64 , αn=-1 )
        p = pars(tpfl,αn)
        (s0,z0,b0)      = p[1]
        (qe,me,ħ)       = p[2]
        (Ms,Msol,tuniv) = p[3]
        (α1,α2)         = p[4]
        (c,yr2s,yr2m)   = p[6]
        (ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI)   = p[5]
        ZER=zero(tpfl)
        if αn < 0
                α = α1
        elseif αn == 0
                α = tpfl(α2)
        else
                α = tpfl(αn)
        end

        if σm>ZER && σet>ZER
            ϑ = (σet/σm)^3
            ξ = ξF(σm,σet,tpfl,αn)
            if ξ<=ZER
                return z0*Ms
            else
                return z0*Ms*σMF(σm,σet,tpfl,αn)
            end
        else
            return z0*Ms
        end
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
"""
    MC1( σm , σet , tpfl=Float64 )
    
This function computes the constraint mass (in solar mass units) 
corresponding to constraint (i) in the article.
"""
function MC1(σm,σet,tpfl=Float64)
    k       = tpfl(1e-15)
    ZER     = zero(tpfl)
    (σmm,σe) = sigtransF(σm,σet,tpfl)
    if σm<ZER || σe<ZER
        return ZER
    else
        return k/σm
    end
end     #---------------------------------------------------------------
