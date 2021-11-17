using DifferentialEquations
using Plots
plotlyjs()

#constants
m = 1.0                             # mass[Kg]
g = 9.81                            # gravitational acceleration [m/s²]

#function to compute the ODEs
function pendulum!(du,u,p,t)
    r,v,θ,ω=u
    du[1]=dr=v                   
    du[2]=dv=r*ω^2+g*cos(θ)
    du[3]=dθ=ω                    
    du[4] =dω=-(g\r)*sin(θ) - (2*v*ω)/r  
end

#initial conditions
r₀=4                            #initial radial distance
v₀=-0.1                          #initial radial velocity 
θ₀ = 1π/18                        # initial angular deflection [rad]
ω₀ = 1                       # initial angular velocity [rad/s]
u₀ = [r₀,v₀,θ₀, ω₀]                   # initial state vector
tspan = (0.0,5.0)             # time interval

#defining the problem and solving
prob = ODEProblem(pendulum!,u₀,tspan)
sol = solve(prob,saveat=0.1)

#plotting
p1=plot(sol,vars=(0,1),label="r",
    color="purple", xaxis="t (s)", yaxis="r (m)")
p2=plot(sol,vars=(0,2),label="v",
    color="blue", xaxis="t (s)", yaxis="v (m/s)")
p3=plot(sol,vars=(0,3),label="θ",
    color="green", xaxis="t (s)", yaxis="θ (rad)")
p4=plot(sol,vars=(0,4),label="ω",
    color="orange", xaxis="t (s)", yaxis="ω (rad/s)")
display(plot(p1,p2,p3,p4,layout=(2,2)))
savefig("./plots/simple length changing pendulum-plot.png")

#trajectory
t=sol.t
U=sol[1:end,:]
(r,v,θ,ω)=[U[x,:] for x in 1:size(U,1)]
x=r.*sin.(θ)
y=zeros(size(t))
z=r.*cos.(θ)
display(plot(x,y,z,linewidth=2,label="trajectory",
    xaxis="x (m)",yaxis="y (m)",zaxis="z (m)",grid=(:on,:black)))
savefig("./plots/simple length changing pendulum-trajectory.png")
