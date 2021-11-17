"""
Julia script for kinematics analysis of four-bar mechanism.
Contains sections for position, velocity, and acceleration analyses.

"""

using Plots
plotlyjs()


# fixed parameters
AB=5
BC=6
CD=5
DA=6
DP=10
θ₁=π

#time limit
t=0:0.01:45

# position analysis

#independent coordinate
ζ₀=π/6 # angle subtended by the A-frame [rad]
β₀=dζ=90*π/(180*45*45)
θ₂=ζ=ζ₀.+β₀.*t.^2

x=-AB.*cos.(θ₂).-DA.*cos.(θ₁)
y=-AB.*sin.(θ₂).-DA.*sin.(θ₁)

z=(x.^2 .+ y.^2 .+ CD^2 .- BC^2)./(2 .* CD)

θ₄=acos.(-z./sqrt.(x.^2 .+ y.^2)) .+ atan.(y,x)  

y_C=CD.*sin.(θ₄)
y_B=AB.*sin.(θ₂)

x_C=CD.*cos.(θ₄).+DA
x_B=AB.*cos.(θ₂)

θ₃=atan.(y_C .- y_B,x_C .- x_B)

#plotting
display(plot(t,θ₂,linewidth=2,xaxis="time (s)",yaxis="θ (rad)",
    label="θ₂",color="blue"))
savefig("./plots/ex3-11.png")
display(plot(t,θ₃,linewidth=2,xaxis="time (s)",yaxis="θ (rad)",
    label="θ₃",color="purple",ylim=[-0.4,0.4]))
savefig("./plots/ex3-12.png")
display(plot(t,θ₄,linewidth=2,xaxis="time (s)",yaxis="θ (rad)",
    label="θ₄",color="orange"))
savefig("./plots/ex3-13.png")

#trajectory of point of suspension
xₚ=DA .+ DP.*cos.(θ₄)
yₚ=DP.*sin.(θ₄)
display(plot(xₚ,yₚ,label="trajectory",xaxis="x (m)",yaxis="y (m)",
    linewidth=2,xlim=[0,15],xticks=0:1:15,ylim=[0,15],yticks=0:1:15))
savefig("./plots/ex3-trajectory.png")


#velocity analysis

#independent coordinate
θ̇₂=dθ₂=2 .*β₀.*t

X=ones(2,length(t))

#solving the equations
for i in 1:length(t)
    A=[BC.*sin.(θ₃[i])  -CD.*sin.(θ₄[i]); BC.*cos.(θ₃[i])  -CD.*cos.(θ₄[i])]
    # X=[θ̇₃;θ̇₄]
    B=[-AB.*θ̇₂[i].*sin.(θ₂[i]); -AB.*θ̇₂[i].*cos.(θ₂[i])]

    Aᴵ=inv(A)
    @assert size(Aᴵ,2)==size(B,1) "dimension mismatch"
    X[:,i]=Aᴵ*B
end

(θ̇₃,θ̇₄)=[X[i,:] for i in 1:size(X,1)]

#plotting
display(plot(t,θ̇₂,linewidth=2,xaxis="time (s)",yaxis="ω (rad/s)",
    label="ω₂",color="blue"))
savefig("./plots/ex3-21.png")
display(plot(t,θ̇₃,linewidth=2,xaxis="time (s)",yaxis="ω (rad/s)",
    label="ω₃",color="purple",ylim=[-0.4,0.4]))
savefig("./plots/ex3-22.png")
display(plot(t,θ̇₄,linewidth=2,xaxis="time (s)",yaxis="ω (rad/s)",
    label="ω₄",color="orange"))
savefig("./plots/ex3-23.png")


#acceleration analysis

#independent coordinate
θ̈₂=dθ̇₂=2*β₀.*ones(length(t))

Y=ones(2,length(t))

#solving the equations
for i in 1:length(t)
    C=[BC*cos(θ₃[i])  -CD*cos(θ₄[i]); BC*sin.(θ₃[i])  -CD*sin.(θ₄[i])]
    # Y=[θ̈₃;θ̈₄]
    D=[(-AB*θ̈₂[i]*cos.(θ₂[i])+AB*(θ̇₂[i]^2)*sin(θ₂[i])+BC*(θ̇₃[i]^2)*sin.(θ₃[i])
            -CD*(θ̇₄[i]^2)*sin(θ₄[i]));
        (-AB*θ̈₂[i]*sin(θ₂[i])-AB*(θ̇₂[i]^2)*cos(θ₂[i])-BC*(θ̇₃[i]^2)*cos.(θ₃[i])
            +CD*(θ̇₄[i]^2)*cos(θ₄[i]))]
    
            Cᴵ=inv(C)
    @assert size(Cᴵ,2)==size(D,1) "dimension mismatch"
    Y[:,i]=Cᴵ*D
end
    
(θ̈₃,θ̈₄)=[Y[i,:] for i in 1:size(Y,1)]

#plotting
display(plot(t,θ̈₂,linewidth=2,xaxis="time (s)",yaxis="α (rad/s^2)",
    label="α₂",color="blue",ylim=[0,0.0025]))
savefig("./plots/ex3-31.png")
display(plot(t,θ̈₃,linewidth=2,xaxis="time (s)",yaxis="α (rad/s^2)",
    label="α₃",color="purple",ylim=[-0.4,0.4]))
savefig("./plots/ex3-32.png")
display(plot(t,θ̈₄,linewidth=2,xaxis="time (s)",yaxis="α (rad/s^2)",
    label="α₄",color="orange",ylim=[0,0.0025]))
savefig("./plots/ex3-33.png")
