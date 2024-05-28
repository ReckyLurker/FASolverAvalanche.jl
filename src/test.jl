using OrdinaryDiffEq, Plots
import ForwardDiff as diff
import LinearAlgebra as LinAlg

# A function describing the first order ODEProblem [second and higher orders need to be decomposed to first order states]
function simplependulum(du, u, p,t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -p[2]*θ - p[1]*dθ
end

# Condition to exit the ODESolver [Required to maintain continuity of the program]
function condition(u,t,integrator)
    sum(abs,integrator(t, Val{1})) - 1e-2
end

# CallbackObject to maintain continuity [no breaks] and exit solver when condition is met. 
cb = ContinuousCallback(condition, terminate!)


# Initial Conditions and System Parameters
d = [0.1,10]
u₀ = [π, 0.5]
tspan = [0.0, 200.0]
prob = ODEProblem(simplependulum, u₀, tspan,d)
sol = solve(prob, Tsit5(), callback = cb)

# Test Function to calculate the settling time for the system
function calculateSettlingTime(d)
    u₀ = [π, 0.5]
    tspan = [0.0, 200.0]
    prob = ODEProblem(simplependulum, u₀, eltype(d).(tspan),d)
    sol = solve(prob, Tsit5(), callback = cb)
    sol.t[end]
end

# Computes the gradient of the settlingtime function.
rd = diff.gradient(calculateSettlingTime, [0.1,10])
r = calculateSettlingTime([0.1, 10])
println(r, " ", rd)
