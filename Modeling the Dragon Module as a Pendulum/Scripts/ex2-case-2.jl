using DifferentialEquations
using Plots
plotlyjs()

function ex2_case_2(IC)

#constants
g=9.81        #acceleration due to gravity [ms^-2]
m=3000        #mass of capsule [Kg]
L=10          #height of A-frame [m]

#initial conditions
r₀,v₀,u₀,traj=IC #traj-trajectory of capsule till now-not needed for calculation
ζ₀=π/6 # angle subtended by the A-frame [rad]
β₀=dζ=90*π/(180*45) 
  
#function to compute lagrangian
function lagrangian!(du,u,p,t)
    θ,ω,ϕ,ψ=u
    r=r₀
    β=dζ=β₀
    ζ=ζ₀+β*t
    du[1]=dθ=ω
    du[2]=dω=(ψ^2)*sin(θ)*cos(θ) + (g*sin(θ)/r) - 
        (L*(β^2)*(cos(ζ)*cos(θ)*cos(ϕ) + sin(ζ)*sin(θ))/r)
    du[3]=dϕ=ψ
    du[4]=dψ=(L*(β^2)*cos(ζ)*sin(ϕ)/(r*sin(θ))) - (2*ω*ψ*cot(θ))
end

#time interval
tspan=(0.0,45.0)

#solving
prob=ODEProblem(lagrangian!,u₀,tspan)
sol=solve(prob,saveat=0.1)

#post-processing
U=sol[1:end,:]
(θ,ω,ϕ,ψ)=[U[x,:] for x in 1:size(U,1)]
t=sol.t
r=5 .* ones(size(t))
β=dζ=β₀
ζ=ζ₀.+β.*t
x=r.*sin.(θ).*cos.(ϕ)-L.*cos.(ζ)
y=r.*sin.(θ).*sin.(ϕ)
z=r.*cos.(θ)+L.*sin.(ζ)
uₜ=U[:,end]
rₜ=r[end]
vₜ=-v₀


# changing coordinate frames (no change)
trajectory=hcat(x,y,z)
FC=(rₜ,vₜ,uₜ,trajectory)

#plotting the solution
 display(plot(x,y,z,label="trajectory",
    xaxis="x (m)",yaxis="y (m)",zaxis="z (m)",
    grid=(:on,:black)))
savefig("./plots/trajectory-2.png")

 p1=plot(sol,vars=(0,1),label="θ",
    color="green", xaxis="t (s)", yaxis="θ (rad)")
 p2=plot(sol,vars=(0,2),label="ω",
    color="orange", xaxis="t (s)", yaxis="ω (rad/s)")
 p3=plot(sol,vars=(0,3),label="ϕ",
    color="purple", xaxis="t (s)", yaxis="ϕ (rad)")
 p4=plot(sol,vars=(0,4),label="ψ",
    color="blue", xaxis="t (s)", yaxis="ψ (rad/s)")
display(plot(p1,p2,p3,p4,layout=(2,2)))
savefig("./plots/plots-2.png")

#energy plot
KE=0.5*m*((L.*β).^2 .+ (r.*ω).^2 .+ (r.*ψ.*sin.(θ)).^2 .+
     2*L.*r.*β.*(ω.*sin.(ζ).*cos.(θ).*cos.(ϕ)
     .-ω.*cos.(ζ).*sin.(θ).-β.*sin.(θ).*sin.(ζ).*sin.(ϕ)))
PE=m*g.*z
TE=KE.+PE
p7=plot(t,KE,label="KE",color="black",
    xaxis="t (s)", yaxis="KE (J)")
p8=plot(t,PE,label="PE",color="green",
    xaxis="t (s)", yaxis="PE (J)")
p9=plot(t,TE,label="TE",color="red",
    xaxis="t (s)", yaxis="TE (J)")
display(plot!(p7,p8,p9))
savefig("./plots/energyplot-2.png")

return FC

end