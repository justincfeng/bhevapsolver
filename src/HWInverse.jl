# Inverse functions for HWFunc.jl
using OrdinaryDiffEq

"""
    dx_dμ(μ, z0)

Compute the derivative dx/dμ where x is defined as:
x = (-z0^(2/3) + (z0^2 + 6μ(9μ + √(3z0^2+81μ^2)))^(1/3))^2/
    (3(z0^(2/3))(z0^2 + 6μ(9μ + √(3z0^2+81μ^2)))^(1/3))
"""
function dx_dμ(μ, z0)
    tpfl = typeof(μ)
    (ONE,TWO,THR,SIX,NIN) = (one(tpfl),tpfl(2),tpfl(3),tpfl(6),tpfl(9))
    ETO = tpfl(81)
    (n1,n2) = (ONE/THR,TWO/THR)
    
    # Helper terms
    sqrt_term = sqrt(THR*z0^2 + ETO*μ^2)
    inner_term = z0^2 + SIX*μ*(NIN*μ + sqrt_term)
    
    # Derivative of sqrt_term
    d_sqrt = (ETO*μ)/(sqrt_term)
    
    # Derivative of inner_term
    d_inner = SIX*(TWO*NIN*μ + sqrt_term + μ*d_sqrt)
    
    # Full derivative using chain rule
    numerator = TWO*(-z0^n2 + inner_term^n1)*inner_term^n1*d_inner
    denominator = THR*(z0^n2)*inner_term^(TWO*n1)
    
    return numerator/denominator
end

"""
    μ_from_τatt(τ_target, z0, μ0=1e-6; abstol=1e-12, reltol=1e-12)

Invert the Δτatt function by solving the ODE dμ/dτ = 1/(dτ/dμ) where
dτ/dμ = (dτ/dx)(dx/dμ) by the chain rule.

Parameters:
- τ_target: The target value of τ to solve for
- z0: The z0 parameter
- μ0: Initial guess for μ (should be small but non-zero)
- abstol: Absolute tolerance for ODE solver
- reltol: Relative tolerance for ODE solver

Returns:
- The value of μ such that Δτatt(μ, z0) ≈ τ_target
"""
function μ_from_τatt(τ_target, z0, μ0=1e-6; abstol=1e-12, reltol=1e-12)
    tpfl = typeof(z0)
    
    # The ODE: dμ/dτ = 1/((dτ/dx)(dx/dμ))
    function dμdτ!(du, u, p, τ)
        μ = u[1]
        # Get x for current μ
        (ONE,TWO,THR,SIX,NIN) = (one(tpfl),tpfl(2),tpfl(3),tpfl(6),tpfl(9))
        ETO = tpfl(81)
        (n1,n2) = (ONE/THR,TWO/THR)
        
        x = (-z0^n2 + (z0^2 + SIX*μ*(NIN*μ + sqrt(THR*z0^2+ETO*μ^2)))^n1)^2/
            (THR*(z0^n2)*(z0^2 + SIX*μ*(NIN*μ +sqrt(THR*z0^2+ETO*μ^2)))^n1)
            
        # Get dτ/dx from dτatt
        dτdx = dτatt(x, z0)
        # Get dx/dμ
        dxdμ = dx_dμ(μ, z0)
        
        # Total derivative using chain rule
        dτdμ = dτdx * dxdμ
        
        # Avoid division by zero
        if abs(dτdμ) < eps(tpfl)
            du[1] = zero(tpfl)
        else
            # dμ/dτ = 1/(dτ/dμ)
            du[1] = one(tpfl)/dτdμ
        end
    end
    
    # Solve from τ=0 to τ_target
    prob = ODEProblem(dμdτ!, [μ0], (tpfl(0), τ_target))
    sol = solve(prob, Vern9(); abstol=abstol, reltol=reltol, save_everystep=false)
    
    return sol.u[end][1]
end

"""
    verify_inverse(μ_test, z0; rtol=1e-6)

Verify the inverse function by comparing μ_from_τatt(Δτatt(μ_test, z0), z0) ≈ μ_test

Parameters:
- μ_test: Test value of μ
- z0: The z0 parameter
- rtol: Relative tolerance for comparison

Returns:
- (is_good, error) where is_good is true if the error is within rtol
"""
function verify_inverse(μ_test, z0; rtol=1e-6)
    τ = Δτatt(μ_test, z0)
    μ_recovered = μ_from_τatt(τ, z0)
    err = abs(μ_recovered - μ_test)/μ_test
    return (err < rtol, err)
end
