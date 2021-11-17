using DifferentialEquations
using Plots
plotlyjs()

# include the case files
include("ex2-case-1.jl")
include("ex2-case-2.jl")
include("ex2-case-3.jl")

#initial conditions
r₀=6
v₀=-0.25
u₀=[17*π/18,0.1,0.1,0.1]

#calling the functions to run them
IC=[r₀,v₀,u₀]
IC2=ex2_case_1(IC)
IC3=ex2_case_2(IC2)
FC=ex2_case_3(IC3);

#overall trajectory
t1=IC2[4]
t2=IC3[4]
t3=FC[4]
t=vcat(t1,t2,t3)
display(plot(t[:,1],t[:,2],t[:,3],label="trajectory",
    xaxis="x (m)",yaxis="y (m)",zaxis="z (m)",grid=(:on,:black)))
savefig("./plots/trajectory.png")

