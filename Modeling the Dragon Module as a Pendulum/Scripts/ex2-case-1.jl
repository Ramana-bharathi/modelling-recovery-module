using DifferentialEquations
using Plots
plotlyjs()

function ex2_case_1(IC)

#constants
g=9.81       #acceleration due to gravity [ms^-2]
m=3000       #mass of capsule [Kg]

#initial conditions
r₀,v₀,u₀=IC
# u₀=[3,0.1,0.1,0.1]
# r₀=6
# v₀=-0.25

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
tspan=(0.0,8.0)

#defining the problem and solving
prob=ODEProblem(lagrangian!,u₀,tspan)
sol=solve(prob,saveat=0.1)

#extracting solution
U=sol[1:end,:]
(θ,ω,ϕ,ψ)=[U[x,:] for x in 1:size(U,1)]
t=sol.t
r=r₀ .+ v₀*t
v=v₀*ones(size(t))
x=r.*sin.(θ).*cos.(ϕ)
y=r.*sin.(θ).*sin.(ϕ)
z=r.*cos.(θ)
uₜ=U[:,end]
rₜ=r[end]
vₜ=v[end]

# # changing coordinate frames
xₚ =x.-(5*(3^0.5))
yₚ=y
zₚ=z.+5
trajectory=hcat(xₚ,yₚ,zₚ)

FC=(rₜ,vₜ,uₜ,trajectory)

#plotting the solution
display(plot(x,y,z,label="trajectory",
    xaxis="x (m)",yaxis="y (m)",zaxis="z (m)",grid=(:on,:black)))
savefig("./plots/trajectory-1.png")

p1=plot(sol.t,r,label="r",
    color="black",xaxis="t (s)", yaxis="r (m)")
p2=plot(sol.t,v,label="v",
    color="red",xaxis="t (s)", yaxis="v (m/s)")
p3=plot(sol,vars=(0,1),label="θ",
    color="green", xaxis="t (s)", yaxis="θ (rad)")
p4=plot(sol,vars=(0,2),label="ω",
    color="orange", xaxis="t (s)", yaxis="ω (rad/s)")
p5=plot(sol,vars=(0,3),label="ϕ",
    color="purple", xaxis="t (s)", yaxis="ϕ (rad)")
p6=plot(sol,vars=(0,4),label="ψ",
    color="blue", xaxis="t (s)", yaxis="ψ (rad/s)")
display(plot(p1,p2,p3,p4,p5,p6,layout=(3,2)))
savefig("./plots/plots-1.png")

#energy plot
KE=0.5*m*(v.^2 .+ (r.*ω).^2 .+ (r.*ψ.*sin.(θ)).^2)
PE=m*g.*z
TE=KE.+PE
p7=plot(t,KE,label="KE",
    color="black",xaxis="t (s)", yaxis="KE (J)")
p8=plot(t,PE,label="PE",
    color="green",xaxis="t (s)", yaxis="PE (J)")
p9=plot(t,TE,label="TE",
    color="red",xaxis="t (s)", yaxis="TE (J)")
display(plot!(p7,p8,p9))
savefig("./plots/energyplot-1.png")

return FC

end