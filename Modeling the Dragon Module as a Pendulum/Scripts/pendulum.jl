using DifferentialEquations
using Plots
plotlyjs()

#constants
l = 1.0                             # length [m]
m = 1.0                             # mass[Kg]
g = 9.81                            # gravitational acceleration [m/s²]

#function to compute ODE
function pendulum!(du,u,p,t)
    du[1] = u[2]                    # θ'(t) = ω(t)
    du[2] = -(g\l)*sin(u[1])  # ω'(t) = -(g/l) sin θ(t)
end

#initial conditions
θ₀ = π/6                         # initial angular deflection [rad]
ω₀ = 0.0                            # initial angular velocity [rad/s]
u₀ = [θ₀, ω₀]                       # initial state vector
tspan = (0.0,100.0)                  # time interval

#defining the problem and solving
prob = ODEProblem(pendulum!,u₀,tspan)
sol = solve(prob,saveat=0.1)

#plotting the solution
display(plot(sol,linewidth=2,xaxis="t (s)",
    color=["green" "blue"],label=["θ [rad]" "ω [rad/s]"],layout=(2,1)))
savefig("./plots/simple pendulum-plot.png")

# trajectory
t=sol.t
U=sol[1:end,:]
(θ,ω)=[U[x,:] for x in 1:size(U,1)]
x=l.*sin.(π.-θ)
y=zeros(size(t))
z=l.*cos.(π.-θ)
display(plot(x,y,z,linewidth=2,label="trajectory",
    xaxis="x (m)",yaxis="y (m)",zaxis="z (m)",grid=(:on,:black)))
savefig("./plots/simple pendulum-trajectory.png")


