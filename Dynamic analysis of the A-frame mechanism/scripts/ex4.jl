"""
Julia script for dynamic analysis of four-bar mechanism.
Includes the kinematic analysis file and computes the reaction
forces and actuating torque required as a function of time.
"""

using Plots
plotlyjs()

include("ex3.jl") #file for kinematic analysis 

K=ex3() #kinematic analysis

#extracting the kinematic state vectors
(θ₂,θ₃,θ₄,θ̇₂,θ̇₃,θ̇₄,θ̈₂,θ̈₃,θ̈₄)=K

#fixed parameters
#lengths and angles
AB=5
BC=6
CD=5
DA=6
DP=10
θ₁=π

#masses
m=3000
m₂=AB*0.2*0.2*7950
m₃=BC*0.15*0.15*7950
m₄=DP*0.25*0.25*7950
g=9.81

#tension force due to the capsule
T=m*g/2

#moment of inertia
I₂=(m₂/12)*AB^2
I₃=(m₃/12)*BC^2
I₄=(m₄/12)*DP^2

#time limit
t=0:0.01:45

#poses and accelerations of centres of gravity of the links
#link 2
x_2=(AB/2).*cos.(θ₂)
ẍ₂=-(AB/2).*θ̈₂.*sin.(θ₂).-(AB/2).*(θ̇₂.^2).*cos.(θ₂)
y_2=(AB/2).*sin.(θ₂)
ÿ₂=(AB/2).*(θ̈₂.*cos.(θ₂).-(θ̇₂.^2).*sin.(θ₂))

#link 3
x_3=AB.*cos.(θ₂).+(BC/2)*cos.(θ₃)
ẍ₃=-AB.*(θ̈₂.*sin.(θ₂).+(θ̇₂.^2).*cos.(θ₂)) .- (BC/2).*(θ̈₃.*sin.(θ₃).+(θ̇₃.^2).*cos.(θ₃))
y_3=AB.*sin.(θ₂).+(BC/2)*sin.(θ₃)
ÿ₃=AB.*(θ̈₂.*cos.(θ₂).-(θ̇₂.^2).*sin.(θ₂)) .+ (BC/2).*(θ̈₃.*cos.(θ₃).-(θ̇₃.^2).*sin.(θ₃))

#link 4
x_4=DA.+(DP/2).*cos.(θ₄)
ẍ₄=-(DP/2).*(θ̈₄.*sin.(θ₄).+(θ̇₄.^2).*cos.(θ₄))
y_4=(DP/2).*sin.(θ₄)
ÿ₄=(DP/2).*(θ̈₄.*cos.(θ₄).-(θ̇₄.^2).*sin.(θ₄))

#locations of the joints
x_A=zeros(size(x_2))
y_A=zeros(size(y_2))
x_B=AB.*cos.(θ₂)
y_B=AB.*sin.(θ₂)
x_C=DA.+CD.*cos.(θ₄)
y_C=CD.*sin.(θ₄)
x_D=DA.*ones(size(x_4))
y_D=y_A=zeros(size(y_4))
x_P=DA.+DP.*cos.(θ₄)
y_P=DP.*sin.(θ₄)

#solution matrix
X=zeros(9,length(t))

#solving the equations
for i in 1:length(t)
    R₁=[1 0 -1 0 0 0 0 0 0]
    R₂=[0 1 0 -1 0 0 0 0 0]
    R₃=[-(y_A[i]-y_2[i]) (x_A[i]-x_2[i]) (y_B[i]-y_2[i]) -(x_B[i]-x_2[i]) 0 0 0 0 1]
    R₄=[0 0 1 0 -1 0 0 0 0]
    R₅=[0 0 0 1 0 -1 0 0 0]
    R₆=[0 0 -(y_B[i]-y_3[i]) (x_B[i]-x_3[i]) (y_C[i]-y_3[i]) -(x_C[i]-x_3[i]) 0 0 0]
    R₇=[0 0 0 0 1 0 1 0 0]
    R₈=[0 0 0 0 0 1 0 1 0]
    R₉=[0 0 0 0 -(y_C[i]-y_4[i]) (x_C[i]-x_4[i]) -(y_D[i]-y_4[i]) (x_D[i]-x_4[i]) 0]
    A=vcat(R₁,R₂,R₃,R₄,R₅,R₆,R₇,R₈,R₉)
    B=[m₂*ẍ₂[i],(m₂*ÿ₂[i]+m₂*g),I₂*θ̈₂[i],m₃*ẍ₃[i],(m₃*ÿ₃[i]+m₃*g),I₃*θ̈₃[i],m₄*ẍ₄[i],(m₄*ÿ₄[i]+m₄*g+T),I₄*θ̈₄[i]+T*(x_P[i]-x_4[i])]
    
    Aᴵ=inv(A)
    @assert size(Aᴵ,2)==size(B,1) "dimension mismatch"
    X[:,i]=Aᴵ*B
end

#extracting the solution vectors
(F₁,F₂,F₃,F₄,F₅,F₆,F₇,F₈,τₐ)=[X[i,:] for i in 1:size(X,1)]

#plotting
display(plot(t,F₁,linewidth=2,xaxis="time (s)",yaxis="F (N)",label="F₁",color="blue"))
savefig("./plots/ex4-41.png")
display(plot(t,F₂,linewidth=2,xaxis="time (s)",yaxis="F (N)",label="F₂",color="green"))
savefig("./plots/ex4-42.png")
display(plot(t,F₃,linewidth=2,xaxis="time (s)",yaxis="F (N)",label="F₃",color="orange"))
savefig("./plots/ex4-43.png")
display(plot(t,F₄,linewidth=2,xaxis="time (s)",yaxis="F (N)",label="F₄",color="purple"))
savefig("./plots/ex4-44.png")
display(plot(t,F₅,linewidth=2,xaxis="time (s)",yaxis="F (N)",label="F₅",color="cyan"))
savefig("./plots/ex4-45.png")
display(plot(t,F₆,linewidth=2,xaxis="time (s)",yaxis="F (N)",label="F₆",color="black"))
savefig("./plots/ex4-46.png")
display(plot(t,F₇,linewidth=2,xaxis="time (s)",yaxis="F (N)",label="F₇",color="brown"))
savefig("./plots/ex4-47.png")
display(plot(t,F₈,linewidth=2,xaxis="time (s)",yaxis="F (N)",label="F₈",color="grey"))
savefig("./plots/ex4-48.png")
display(plot(t,τₐ,linewidth=2,xaxis="time (s)",yaxis="M (Nm)",label="τₐ",color="red"))
savefig("./plots/ex4-49.png")
