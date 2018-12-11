#-----------------------------------------------------------#
# Nonlinear ACOPF problem
# Last updated: 1/30/2016
# Created by Harsha Nagarajan at LANL
#-----------------------------------------------------------#


# Define packages to be used:
using JuMP, Distributions, MathProgBase, Ipopt, JuMPChance
using Gurobi



include("acopf_input.jl")
include("model_nl.jl")
include("jacobians.jl")
include("model_nl_penalty.jl")

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

#line_limits[119]=10*line_limits[119];
#println("Limit on line 119:  ",line_limits[119])


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

println("bus: ", numbuses)
println("line: ",numlines)
println("gen: ", numgenerators)


# number of iterations for sampling
N_iter = 2


total_violation_mean = []
total_violation_nonzero =[]
expected_cost = []
std_cost = []

delta_pmax_iter =[]


# Solve the CC-ACOPF for different values of parameter ϵ
for ϵ in [0.2 ] #0.05 0.01 0.005 0.001 0.0005 0.0001]

status_det, cost_det, status_cc, cost_cc, v_cc, gp_cc, gq_cc, α_cc, γ_cc, r_cc,  r_det, γ_det, total_dev,
 	gp_det, gq_det, v_det, θ_det, time_det  =
 	createandsolvemodel(thermalLimitscale, generators, buses, lines, farm, ϵ,voll)

	println("Det: ", status_det, "  Cost: ", round(cost_det*10)/10, " Time (s):", time_det)
	println("CC:  epsilon: ", ϵ  ," ", status_cc, "  Cost: ", round(cost_cc*10)/10)

end


#	println("Number of generators in the solution: ", length(gp_det))
#	println("Number of generators in the input data: ", length(generators))



srand(100)

objective_iter=[]
status_iter=[]
s_validation_iter =[]
vmin_validation_iter =[]
vmax_validation_iter =[]
cost_iter =[]

# Run a feasibility test for a given value of epsilon

for iter in 1:N_iter

println("Value of epsilon: ", ϵ, "   Value of iteration: ", iter)

farm_error = []
farm_p_real = []
farm_q_real = []
slack_pmax = zeros(1,numgenerators)

# Generate the real-time output of wind generators using random sampling
for f in 1:length(farm)
	distribution = Normal(0.0, farm[f].σ)
	sample = rand(distribution)
	push!(farm_error, sample)
	push!(farm_p_real, farm[f].μ + sample)
	push!(farm_q_real, γ_cc[f]*farm_p_real[end])
end

total_error = sum(farm_error) # omega in the paper as used in eq. 7 (same convention on alpha)
println("Total deviation: ", total_error)
gp_real =  gp_cc - α_cc * total_error # the value of gp_real is meaningless for the slack bus and not used in acopf_for_validation below
# gp_real =  gp_cc - 1/numgenerators* total_error

# compute the generation slack
for i in 1:numgenerators
	slack_pmax[i] =  generators[i].Pgmax -gp_real[i] + 0
end

gq_real =  gq_cc
v_real = v_cc

push!(delta_pmax_iter, mean(slack_pmax))


cost_real = sum(generators[i].pi1*gp_real[i]^2 + generators[i].pi2*gp_real[i] + generators[i].pi3 for i=1:numgenerators)
println("Show the cost unique:", cost_real)
push!(cost_iter, cost_real)

# Recompute the power flow to assesses the solution
status_validation, objective_validation,  slack_s_validation, slack_vmin_validation, slack_vmax_validation =
	acopf_for_validation(thermalLimitscale, generators, buses, lines, farm, gp_real, gq_real, v_real, farm_p_real, farm_q_real, 100, 100,100)

push!(objective_iter, objective_validation)
push!(status_iter, status_validation)
push!(s_validation_iter, slack_s_validation)
push!(vmin_validation_iter, slack_vmin_validation)
push!(vmax_validation_iter, slack_vmax_validation)



end

push!(expected_cost, mean(cost_iter))
#push!(std_cost, std(cost_iter))
push!(total_violation_mean, mean(objective_iter))
push!(total_violation_nonzero, sum(objective_iter .>= 1e-3) + sum(status_iter .!= :Optimal))


	s_violation_count = zeros(1, numlines)
	vmin_violation_count = zeros(1, numbuses)
	vmax_violation_count = zeros(1, numbuses)


	println(s_violation_count)

	for iter in 1:N_iter
		println("iter: ", iter, "  Status: ", status_iter[iter], "  Objective: ", objective_iter[iter])

		for l in 1:numlines
			if s_validation_iter[iter][l]>1e-3
				println("Iter: ", iter, "  Line: ", l, " Violation: ", s_validation_iter[iter][l])
				s_violation_count[l] = s_violation_count[l] +1
			end
		end

		for b in 1:numbuses
			if vmin_validation_iter[iter][b]>1e-3
				println("Iter: ", iter, "  Bus ", b, " Vmin violation: ", vmin_validation_iter[iter][b])
				vmin_violation_count[b] = vmin_violation_count[b] + 1
			end


			if vmax_validation_iter[iter][b]>1e-3
				println("Iter: ", iter, "  Bus ", b, " Vmax violation: ", vmax_validation_iter[iter][b])
				vmax_violation_count[b] = vmax_violation_count[b] + 1
			end

		end
	end





	println("******************************")
	println("Print violations per each line")
	println("******************************")
	for l in 1:numlines
		if s_violation_count[l]>0
		println("Epsilon: ",  ϵ, " Line: ", l, " Total violation: ", s_violation_count[l])
		end
	end

	println("******************************")
	println("Print violations per each bus")
	println("******************************")


	for b in 1:numbuses
		if vmin_violation_count[b]>0
		println("Epsilon: ",  ϵ, " Bus: ", b, " Vmin total violation: ", vmin_violation_count[b])
		end

		if vmax_violation_count[b]>0
		println("Epsilon: ",  ϵ, " Bus: ", b, " Vmax total violation: ", vmax_violation_count[b])
		end

	end



end
