 
#-----------------------------------------------------------------------
module HWslv      # HWs module
#-----------------------------------------------------------------------

using LinearAlgebra

#-----------------------------------------------------------------------
#   NUMBERS
#-----------------------------------------------------------------------
"""
    nums( tpfl=Float64 )
    
This function recomputes numbers in the precision specified in the 
argument.
"""
function nums( tpfl )
    (ONE,TWO,THR,FOU,FIV) = (tpfl(1),tpfl(2),tpfl(3),tpfl(4),tpfl(5))
    (SIX,SEV,EIG,NIN,PI)  = (tpfl(6),tpfl(7),tpfl(8),tpfl(9),tpfl(π))
    return (ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   muattr
#-----------------------------------------------------------------------
"""
    muattr( Y , p )
    
This function computes the approximate attractor curve, written as a 
function `μ(Y)`. The inputs are `Y` and the parameter tuple `p`, the
latter being the output of the `paramset` function.
"""
function muattr( Y , p )
    (s0,z0,b0)   = p
    tpfl  = typeof(Y)
    ONE   = tpfl(1)

    return (√(Y)*z0)/((ONE + √(tpfl(ONE-Y)))^2)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   Yattr
#-----------------------------------------------------------------------
"""
    Yattr( μ , p )
    
This function computes the approximate attractor curve, written as a 
function `Y(μ)`. The inputs are `μ` and the parameter tuple `p`, the
latter being the output of the `paramset` function.
"""
function Yattr( μ , p )
    (s0,z0,b0)   = p
    tpfl  = typeof(Y)
    (ONE,THR,FOU,SIX,NIN) = (tpfl(1),tpfl(3),tpfl(4),tpfl(6),tpfl(9))
    TTS = tpfl(27)
    z0s = z0^2
    h = (z0s + SIX*μ*(NIN*μ + sqrt(THR)*sqrt(z0s + TTS*μ^2)))^(1/3)
    x = (h - z0s^(1/3))^2/(THR*h*z0s^(1/3))
    return FOU*x/((ONE+x)^2)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   YH
#-----------------------------------------------------------------------
"""
    YH( μ , Y1 )

This function provides an approximate phase space solution for the 
system in the Hawking phase.
"""
function YH( μ , Y1 )
    return Y1/μ^2
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   YS
#-----------------------------------------------------------------------
"""
    YS( μ , μ0 )

This function provides an approximate phase space solution for the 
system in the Schwinger phase.
"""
function YS( μ , μ0 )
    return ((2*μ-μ0)*μ0)/μ^2
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   SFunc
#-----------------------------------------------------------------------
"""
    SFunc( X , p )
    
This function evaluates the exponential factor in the simplified Hiscock
and Weems equations, given the input vector `X=[μ,Y]` (where `μ:=M/Ms` 
and `Y:=Q/M`) and parameter tuple `p`, the latter being the output of
the `paramset` function.
"""
function SFunc( X , p )
    μ = X[1]
    Y = X[2]
    (s0,z0,b0)   = p
    tpfl  = typeof(Y)
    ONE   = tpfl(1)
 
    return exp( b0*tpfl(z0-μ*((ONE + √(tpfl(ONE-Y)))^2)/√Y) )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   HFunc
#-----------------------------------------------------------------------
"""
    HFunc( X , p )
    
This function evaluates the Hawking factor in the simplified Hiscock and
Weems equations, given the input vector `X=[μ,Y]` (where `μ:=M/Ms` 
and `Y:=Q/M`) and parameter tuple `p`, the latter being the output of the
`paramset` function.
"""
function HFunc( X , p )
    μ = X[1]
    Y = X[2]
    (s0,z0,b0)   = p
    tpfl  = typeof(Y)
    (ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI) = nums(tpfl)

    Ξ = √(tpfl(ONE-Y))
    Θ = √(tpfl(NIN-EIG*Y))

    return ((THR + Θ)^4*(ONE - Y)^2)/((ONE + Ξ)^4*(THR-TWO*Y+Θ)*μ^2)
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   F
#-----------------------------------------------------------------------
"""
    F( X , p )
    
This function computes `[ μ̇ , Ẏ ]` for the Hiscock and Weems equation
(where `μ:=M/Ms` and `Y:=Q/M`), given the input vector `X=[μ,Y]` and
parameter tuple `p`, the latter being the output of the `paramset`
function.
"""
function F( X , p , SF=SFunc , HF=HFunc )
    μ = X[1]
    Y = X[2]
    (s0,z0,b0)   = p
    tpfl  = typeof(Y)
    (ONE,TWO,THR,FOU,FIV,SIX,SEV,EIG,NIN,PI) = nums(tpfl)

    Ξ = √(tpfl(ONE-Y))
 
    S    = SF( X , p )
    H    = HF( X , p )
    
    dμ    = -(H+S*Y^2)/(ONE+Ξ)^4
    dY    = TWO*Y*(H - S*(1-Y+Ξ)*Y) /(μ*(ONE+Ξ)^4)

    return [ dμ , dY ]
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
#   Frev
#-----------------------------------------------------------------------
"""
    Frev( X , p )
    
This function is the same as `F` except for a sign. This can be used to 
evolve the system backwards in time.
"""
function Frev( X , p )
    return -F( X , p )
end     #---------------------------------------------------------------

#-----------------------------------------------------------------------
end     # End scope of module HWslv
#-----------------------------------------------------------------------