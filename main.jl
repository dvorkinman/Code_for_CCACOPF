#-----------------------------------------------------------#
# Nonlinear ACOPF problem
# Last updated: 1/30/2016
# Created by Harsha Nagarajan at LANL
#-----------------------------------------------------------#


# Define packages to be used:
using JuMP, Distributions, MathProgBase, Ipopt, JuMPChance

ϵ = 0.01

include("acopf_input.jl")
include("model_nl.jl")
include("jacobians.jl")


# Read the input data:
path = pwd()"/case.dat"
casefilename, refbus, loadscale, thermalLimitscale, extras, mvaBase = readconfig(path)

# Read specific data for the test case:
generators, generatorlist, buses, lines = readcase(casefilename, extras, loadscale, thermalLimitscale, mvaBase)

farm = Farm[]

numbuses = length(buses)
numlines = length(lines)
numgenerators = length(generators)

# If needed, print out the bus and gen data
#	println(">>>> Buses: $(numbuses),  Lines: $(numlines)")
# 	println(">>>> Generators: $(numgenerators)")
# 	for i in 1:numgenerators
# 		println("Number", i,"P max in MW: ",  generators[i].Pgmax*100)
# 	end

line_limits= [ 175	175	500	175	175	175	500	500	500	175	175	175	175	175	175	175	175	175	175	175	500	175	175	175	175	175	175	175	175	175	500	500	500	175	175	500	175	500	175	175	140	175	175	175	175	175	175	175	175	500	500	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	500	175	175	500	500	500	500	500	500	500	175	175	500	175	500	175	175	500	500	175	175	175	175	175	175	175	500	175	175	175	175	175	175	500	500	175	500	500	200	200	175	175	175	500	500	175	175	500	500	500	175	500	500	175	175	175	175	175	175	175	175	175	175	200	175	175	175	175	175	175	175	175	175	500	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	500	175	175	175	500	175	175	175]
for i in 1:length(lines)
	lines[i].u = 0.99*thermalLimitscale * line_limits[i]/mvaBase
end

wp = 1.25
factor_σ =  1.25*wp
voll = 10000



push!(farm, Farm(70.0/100*wp, 		factor_σ *7.0/100, 	3))
push!(farm, Farm(147.0/100*wp, 	factor_σ *14.7/100, 	8))
push!(farm, Farm(102.0/100*wp, 	factor_σ *10.2/100, 	11))
push!(farm, Farm(105.0/100*wp, 	factor_σ *10.5/100, 	20))
push!(farm, Farm(113.0/100*wp, 	factor_σ *11.3/100, 	24))
push!(farm, Farm(84.0/100*wp,      factor_σ *8.4/100, 	26))
push!(farm, Farm(59.0/100*wp,  	factor_σ * 5.9/100, 	31))
push!(farm, Farm(250.0/100*wp, 	factor_σ *25.0/100, 	38))
push!(farm, Farm(118.0/100*wp,  	factor_σ *11.8/100, 	43))
push!(farm, Farm(76.0/100*wp,  	factor_σ *7.6/100, 	49))
push!(farm, Farm(72.0/100*wp,  	factor_σ *7.2/100, 	53))

status_det, cost_det, status_cc, cost_cc,  gp_cc, gq_cc, α_cc, γ_cc, r_cc,  gp_det, r_det, γ_det, total_dev,
 	gp_initial, gq_initial, v_initial, θ_initial, time_det  =
 	createandsolvemodel(thermalLimitscale, generators, buses, lines, farm, ϵ,voll)

# println("Det: ", status_det, "  Cost: ", round(cost_det*10)/10, " Time (s):", time_det)
# println("CC:  ", status_cc, "  Cost: ", round(cost_cc*10)/10)
