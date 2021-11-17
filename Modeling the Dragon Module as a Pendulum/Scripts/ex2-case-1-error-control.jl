using DifferentialEquations
using Plots
plotlyjs()



#constants
g=9.81       #acceleration due to gravity [ms^-2]
m=3000       #mass of capsule [Kg]

#initial conditions
u₀=[3,0.1,0.1,0.5]
r₀=6
v₀=0

#function to compute lagrangian
function lagrangian!(du,u,p,t)
    θ,ω,ϕ,ψ=u
    v=v₀
    r=r₀+v*t
    du[1]=dθ=ω
    du[2]=dω=sin(θ)*cos(θ)*ψ^2-((2*v*ω)/r)+((g*sin(θ))/r)
    du[3]=dϕ=ψ
    du[4]=dψ=-((2*v*ψ)/r)-(2*cot(θ)*ω*ψ)
end

#time interval
tspan=(0.0,20.0)

#defining the problem and solving
prob=ODEProblem(lagrangian!,u₀,tspan)
sol1=solve(prob,saveat=0.1,reltol=1e-1,abstol=1e-4)
sol2=solve(prob,saveat=0.1)
sol3=solve(prob,saveat=0.1,reltol=1e-5,abstol=1e-8)

#extracting solution
U1=sol1[1:end,:]
(θ1,ω1,ϕ1,ψ1)=[U1[x,:] for x in 1:size(U1,1)]
t=sol1.t

U2=sol2[1:end,:]
(θ2,ω2,ϕ2,ψ2)=[U2[x,:] for x in 1:size(U2,1)]
t=sol2.t

U3=sol3[1:end,:]
(θ3,ω3,ϕ3,ψ3)=[U3[x,:] for x in 1:size(U3,1)]
t=sol3.t


#plotting the solutions for different tolerances

p1=plot(sol1,vars=(0,1),label="rtol=1e-1,atol=1e-4",
    color="green", xaxis="t (s)", yaxis="θ (rad)")
p2=plot!(p1,sol2,vars=(0,1),label="rtol=1e-3,atol=1e-6",
    color="blue", xaxis="t (s)", yaxis="θ (rad)")
p3=plot!(p2,sol3,vars=(0,1),label="rtol=1e-5,atol=1e-8",
    color="red", xaxis="t (s)", yaxis="θ (rad)")
display(plot(p3))
savefig("./plots/error_control_1.png")

p4=plot(sol1,vars=(0,3),label="rtol=1e-1,atol=1e-4",
    color="green", xaxis="t (s)", yaxis="ϕ (rad)")
p5=plot!(p4,sol2,vars=(0,3),label="rtol=1e-3,atol=1e-6",
    color="blue", xaxis="t (s)", yaxis="ϕ (rad)")
p6=plot!(p5,sol3,vars=(0,3),label="rtol=1e-5,atol=1e-8",
    color="red", xaxis="t (s)", yaxis="ϕ (rad)")
display(plot(p6))
savefig("./plots/error_control_2.png")
