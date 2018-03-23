function createandsolvemodel(thermalLimitscale, generators, buses, lines, farms, ϵ, voll)

    numbuses = length(buses)
    numlines = length(lines)
    numgenerators = length(generators)
    numfarms = length(farms)
    θᵘ = deg2rad(float(extras["theta_u"]))

    m = Model(solver=IpoptSolver())


    PV = find(map(b->b.kind==:PV, buses))
    PQ = find(map(b->b.kind==:PQ, buses))
    REF = find(map(b->b.kind==:Ref, buses))
    @assert length(REF) == 1

    total_var = sum([f.σ^2 for f in farms])
    total_deviation = quantile(Normal(0,sqrt(total_var)),1-ϵ)
    @show total_deviation

    # Continuous Variables
    @variable(m, pij[1:numlines])
    @variable(m, qij[1:numlines])
    @variable(m, pji[1:numlines])
    @variable(m, qji[1:numlines])
    @variable(m, θ[1:numbuses])
    @variable(m, -θᵘ <= θij[1:numlines] <= θᵘ)
    @variable(m, gp[1:numgenerators])
    @variable(m, gq[1:numgenerators])
    @variable(m, r[1:numgenerators] >= 0)
    @variable(m, v[1:numbuses] >= 0)
    @variable(m, γ[1:numfarms] >= 0, start=0.5)


    @constraint(m, sum(r) >= total_deviation)

    # Objective
    @NLobjective(m, Min, sum{generators[i].pi1*gp[i]^2 + generators[i].pi2*gp[i] + generators[i].pi3, i=1:numgenerators})

    # AC equations
    pij_constr = []
    qij_constr = []
    pji_constr = []
    qji_constr = []

    for i=1:numlines
        @constraint(m, θij[i] == (θ[lines[i].head] - θ[lines[i].tail]))

        # flow leaving i towards j
        @NLconstraint(m, con, pij[i] == (1/lines[i].ratio)^2 * lines[i].γ * v[lines[i].head]^2 - (1/lines[i].ratio)*lines[i].γ * v[lines[i].head]*v[lines[i].tail] *cos(θ[lines[i].head] - θ[lines[i].tail]) - (1/lines[i].ratio)*lines[i].β * v[lines[i].head]*v[lines[i].tail] *sin(θ[lines[i].head] - θ[lines[i].tail]))
        push!(pij_constr, con)

        @NLconstraint(m, con, qij[i] == -(1/lines[i].ratio)^2*(lines[i].β + lines[i].b_charge/2)* v[lines[i].head]^2 + (1/lines[i].ratio)*lines[i].β * v[lines[i].head]*v[lines[i].tail] *cos(θ[lines[i].head] - θ[lines[i].tail]) - (1/lines[i].ratio)*lines[i].γ * v[lines[i].head]*v[lines[i].tail] *sin(θ[lines[i].head] - θ[lines[i].tail]))
        push!(qij_constr, con)

        # flow leaving j towards i
        @NLconstraint(m, con, pji[i] ==  lines[i].γ * v[lines[i].tail]^2 - (1/lines[i].ratio) * lines[i].γ * v[lines[i].head]*v[lines[i].tail] *cos(θ[lines[i].head] - θ[lines[i].tail]) + (1/lines[i].ratio) * lines[i].β * v[lines[i].head]*v[lines[i].tail] *sin(θ[lines[i].head] - θ[lines[i].tail]))
        push!(pji_constr, con)

        @NLconstraint(m, con, qji[i] == -(lines[i].β + lines[i].b_charge/2) * v[lines[i].tail]^2 +(1/lines[i].ratio) * lines[i].β * v[lines[i].head]*v[lines[i].tail] *cos(θ[lines[i].head] - θ[lines[i].tail]) + (1/lines[i].ratio) * lines[i].γ * v[lines[i].head]*v[lines[i].tail] * sin(θ[lines[i].head] - θ[lines[i].tail]))
        push!(qji_constr, con)

        # Line limits
        if lines[i].u > 0
            @NLconstraint(m, pij[i]^2 + qij[i]^2 <= lines[i].u^2)
            @NLconstraint(m, pji[i]^2 + qji[i]^2 <= lines[i].u^2)
        end
    end

    # Flow balance

    @expression(m, sum_pgen[i=1:numbuses], sum{gp[k], k in buses[i].genids} + sum{farms[k].μ, k=1:numfarms; farms[k].bus==i})
    @expression(m, sum_qgen[i=1:numbuses], sum{gq[k], k in buses[i].genids} + sum{γ[k]*farms[k].μ, k=1:numfarms; farms[k].bus==i})
    for i=1:numbuses
        @constraint(m, sum_pgen[i] - buses[i].Pd - buses[i].Gs * v[i]^2 == sum{pij[k], k in buses[i].outlist} + sum{pji[k], k in buses[i].inlist})
        @constraint(m, sum_qgen[i] - buses[i].Qd + buses[i].Bs * v[i]^2== sum{qij[k], k in buses[i].outlist} + sum{qji[k], k in buses[i].inlist})
        if buses[i].kind == :Ref
            @constraint(m, θ[i] == 0)
        elseif buses[i].kind == :PV
        end
    end

    # Operational limits with reserves
    @constraint(m, c_gen1[i=1:numgenerators], gp[i] - r[i] >= generators[i].Pgmin)
    @constraint(m, c_gen2[i=1:numgenerators], gp[i] + r[i] <= generators[i].Pgmax)
    @constraint(m, c_gen3[i=1:numgenerators], gq[i] >= generators[i].Qgmin)
    @constraint(m, c_gen4[i=1:numgenerators], gq[i] <= generators[i].Qgmax)

    @constraint(m, c_v1[i=1:numbuses], v[i] <= buses[i].Vmax)
    @constraint(m, c_v2[i=1:numbuses], v[i] >= buses[i].Vmin)

    tic()
    status_det = solve(m)

    time_det = [] # getsolvetime(m::Model)

    gp_initial = [] #Array(Float64, numgenerators, 1)
    gq_initial = [] #Array(Float64, numgenerators, 1)
    v_initial = [] #Array(Float64, numbuses, 1)
    θ_initial = [] #Array(Float64, numbuses, 1)


    gp_initial = getvalue(gp)
    gq_initial = getvalue(gq)
    v_initial = getvalue(v)
    θ_initial = getvalue(θ)
    # f_opt = getvalue(pij)^2 + getvalue(qij)^2

    println("*******************************************")
    println("Objective value: ", getobjectivevalue(m))
    println("Time taken (seconds): ", toq())
    println("Phase angle bound (deg-ree): ", rad2deg(θᵘ))
    println("*******************************************")

    # linearized map from (v,θ) to flows
    x = [v;θ]
    flows = [pij;pji;qij;qji]
    flows_constr = [pij_constr;pji_constr;qij_constr;qji_constr]

    vθ_to_flows_const, vθ_to_flows_J = extract_jacobian(m, flows_constr, x)
    vθ_to_flows_J *= -1 # JuMP moves lhs to rhs
    vθ_to_flows_const -= getvalue(flows)
    vθ_to_flows_const *= -1
    vθ_to_flows_const -= vθ_to_flows_J*getvalue(x)

    # flows ≈ vθ_to_flows_const + vθ_to_flows_J*(x)
    # @assert maximum(abs(getvalue(flows) - (vθ_to_flows_const+vθ_to_flows_J*getvalue(x)))) <= 1e-6

    # now set up map from (v,θ) to nodal (p,q)
    vθ_to_node_p_J = zeros(numbuses,2*numbuses)
    vθ_to_node_p_const = zeros(numbuses)
    vθ_to_node_q_J = zeros(numbuses,2*numbuses)
    vθ_to_node_q_const = zeros(numbuses)
    for i=1:numbuses
        # p = Pd + Gs*v[i]^2 + sum{pij[k], k in buses[i].outlist} + sum{pji[k], k in buses[i].inlist}
        vθ_to_node_p_const[i] = buses[i].Pd
        vθ_to_node_p_const[i] -= buses[i].Gs*getvalue(v[i])^2
        vθ_to_node_p_J[i,i] += 2*buses[i].Gs*getvalue(v[i])
        for k in buses[i].outlist
            vθ_to_node_p_J[i,:] += vθ_to_flows_J[k,:] # pij
            vθ_to_node_p_const[i] += vθ_to_flows_const[k]
        end
        for k in buses[i].inlist
            vθ_to_node_p_J[i,:] += vθ_to_flows_J[k+numlines,:] # pji
            vθ_to_node_p_const[i] += vθ_to_flows_const[k+numlines]
        end

        # q = buses[i].Qd - buses[i].Bs * v[i]^2 + sum{qij[k], k in buses[i].outlist} + sum{qji[k], k in buses[i].inlist})
        vθ_to_node_q_const[i] = buses[i].Qd
        vθ_to_node_q_const[i] += buses[i].Bs*getvalue(v[i])^2
        vθ_to_node_q_J[i,i] -= 2*buses[i].Bs*getvalue(v[i])
        for k in buses[i].outlist
            vθ_to_node_q_J[i,:] += vθ_to_flows_J[k+2*numlines,:] # qij
            vθ_to_node_q_const[i] += vθ_to_flows_const[k+2*numlines]
        end
        for k in buses[i].inlist
            vθ_to_node_q_J[i,:] += vθ_to_flows_J[k+3numlines,:] # qji
            vθ_to_node_q_const[i] += vθ_to_flows_const[k+3numlines]
        end
    end

    # p ≈ vθ_to_node_p_const + vθ_to_node_p_J*(x)
    # q ≈ vθ_to_node_q_const + vθ_to_node_q_J*(x)

    pgen = [sum(Float64[getvalue(gp[k]) for k in buses[i].genids]) for i in 1:numbuses]
    qgen = [sum(Float64[getvalue(gq[k]) for k in buses[i].genids]) for i in 1:numbuses]
    for k in 1:numfarms
        pgen[farms[k].bus] += farms[k].μ
        qgen[farms[k].bus] += getvalue(γ[k])*farms[k].μ
    end

    @show maximum(pgen - (vθ_to_node_p_const + vθ_to_node_p_J*getvalue(x)))
    @show maximum(qgen - (vθ_to_node_q_const + vθ_to_node_q_J*getvalue(x)))

    vθ_to_node_pq_J = [ vθ_to_node_p_J
                        vθ_to_node_q_J ]
    vθ_to_node_pq_const = [ vθ_to_node_p_const
                            vθ_to_node_q_const ]

    #=
    [ y ] = [ A B ][ η ] + [ c₁ ]
    [ z ]   [ C D ][ β ]   [ c₂ ]
    where η = [ v_pq, θ_pq, θ_pv ]
          β = [ v_pv, v_θv, θ_θv ]
          y = [ p_pq, p_pv, q_pq ]
          z = [ p_θv, q_pv, q_θv ].
    Want to express η as a function of y, β,
        so η = A^{-1}(y - B*β - c₁)
    Then z = C*η + D*β + c₂ is easy to get as well.
    =#

    y_rows = [PQ;PV;PQ+numbuses]
    z_rows = [REF;PV+numbuses;REF+numbuses]
    η_cols = [PQ;PQ+numbuses;PV+numbuses]
    β_cols = [PV;REF;REF+numbuses]

    A = sparse(vθ_to_node_pq_J[y_rows,η_cols])
    #B = sparse(vθ_to_node_pq_J[y_rows,β_cols])
    C = sparse(vθ_to_node_pq_J[z_rows,η_cols])
    #D = sparse(vθ_to_node_pq_J[z_rows,β_cols])
    #c₁ = vθ_to_node_pq_const[y_rows]
    #c₂ = vθ_to_node_pq_const[z_rows]
    Ainv = A\eye(size(A,1))
    Ainv[abs(Ainv) .<= 1e-6] = 0.0

    Ainvsp = sparse(Ainv)
    @show nnz(Ainvsp)
    @show size(A,1)^2

    m_chance = ChanceModel(solver=IpoptSolver())

    ## Deterministic part of the model

    @variable(m_chance, θ_chance[1:numbuses]) # nominal angles
    @variable(m_chance, gp_chance[1:numgenerators])
    @variable(m_chance, generators[i].Qgmin <= gq_chance[i=1:numgenerators] <= generators[i].Qgmax)
    @variable(m_chance, v_chance[i=1:numbuses]) # nominal voltage
    @variable(m_chance, r_chance[1:numgenerators] >= 0)
    @variable(m_chance, pij_chance[1:numlines])
    @variable(m_chance, qij_chance[1:numlines])
    @variable(m_chance, pji_chance[1:numlines])
    @variable(m_chance, qji_chance[1:numlines])
    @variable(m_chance, γ_chance[1:numfarms] >= 0)

    @variable(m_chance, genquad[i=1:numgenerators] >= 0)
    @variable(m_chance, genquad_aux[i=1:numgenerators] == 1)
    @constraint(m_chance, genquad_rsoc[i=1:numgenerators],
                    genquad[i]*genquad_aux[i] >= gp_chance[i]^2)
    # Objective
    @objective(m_chance, Min, sum{generators[i].pi1*genquad[i] + generators[i].pi2*gp_chance[i] + generators[i].pi3, i=1:numgenerators})

    # add linearization constraints for the deterministic part
    @expression(m_chance, sum_pgen_chance[i=1:numbuses], sum{gp_chance[k], k in buses[i].genids} + sum{farms[k].μ, k=1:numfarms; farms[k].bus==i})
    @expression(m_chance, sum_qgen_chance[i=1:numbuses], sum{gq_chance[k], k in buses[i].genids} + sum{γ_chance[k]*farms[k].μ, k=1:numfarms; farms[k].bus==i})
    # p ≈ vθ_to_node_p_const + vθ_to_node_p_J*(x)
    # q ≈ vθ_to_node_q_const + vθ_to_node_q_J*(x)
    x_chance = [v_chance;θ_chance]
    @constraint(m_chance, sum_pgen_chance .== vθ_to_node_p_const + vθ_to_node_p_J*x_chance)
    @constraint(m_chance, sum_qgen_chance .== vθ_to_node_q_const + vθ_to_node_q_J*x_chance)

    @constraint(m_chance, sum(r_chance) >= total_deviation)

    # Operational limits with reserves
    @constraint(m_chance, [i=1:numgenerators], gp_chance[i] - r_chance[i] >= generators[i].Pgmin)
    @constraint(m_chance, [i=1:numgenerators], gp_chance[i] + r_chance[i] <= generators[i].Pgmax)

    # [pij,pji,qij,qji] ≈ vθ_to_flows_const + vθ_to_flows_J*(x)
    @constraint(m_chance, [pij_chance;pji_chance;qij_chance;qji_chance] .==
                            vθ_to_flows_const + vθ_to_flows_J*x_chance)
    # deterministic line limits
    for i in 1:numlines
        if lines[i].u > 0
            @constraint(m_chance, pij_chance[i]^2 + qij_chance[i]^2 <= lines[i].u^2)
            @constraint(m_chance, pji_chance[i]^2 + qji_chance[i]^2 <= lines[i].u^2)
        end
    end

    @constraint(m_chance, θ_chance[REF] .== 0)

    ## Response part of the model

    @variable(m_chance, α_chance[1:numgenerators] >= 0)
    @constraint(m_chance, sum(α_chance) == 1)

    @indepnormal(m_chance, ω[i=1:numfarms], mean=0, var=farms[i].σ^2)
    Ω = sum(ω)
    @expression(m_chance, α_by_bus[b=1:numbuses], sum{ α_chance[i], i in buses[b].genids })
    haswind = falses(numbuses)
    for k in 1:numfarms
        haswind[farms[k].bus] = true
    end

    p_response = Vector{JuMPChance.CCAffExpr}(numbuses)
    q_response = Vector{JuMPChance.CCAffExpr}(numbuses)
    for i in 1:numbuses
        @expression(m_chance, sum_presponse, sum{ω[k], k=1:numfarms; farms[k].bus==i} - sum{α_by_bus[i]*Ω; haswind[i]})
        @expression(m_chance, sum_qresponse, sum{γ_chance[k]*ω[k], k=1:numfarms; farms[k].bus==i})
        p_response[i] = AffExpr() + sum_presponse
        q_response[i] = AffExpr() + sum_qresponse
    end

    # first solve for Δη = Δ[v_pq,θ_pq,θ_pv] given Δy = Δ[p_pq, p_pv, q_pq] and Δβ = Δ[v_pv,v_θv,θ_θv].
    # Actually Δβ is zero since they're fixed, don't respond.
    # so Δη = A^{-1}(Δy)

    Δy = JuMPChance.CCAffExpr[ p_response[PQ]; p_response[PV]; q_response[PQ] ]
    tic()
    Δη = Ainvsp*Δy

    v_response = zeros(JuMPChance.CCAffExpr,numbuses) + v_chance
    v_response[PQ] += Δη[1:length(PQ)]
    θ_response = zeros(JuMPChance.CCAffExpr,numbuses) + θ_chance
    θ_response[PQ] += Δη[(length(PV)+1):(length(PV)+length(PQ))]
    θ_response[PV] += Δη[(length(PV)+length(PQ)+1):(2length(PV)+length(PQ))]
    toc()
    # impose chance constraints on v and θ
    # some of these are deterministic
    v_limits = @constraint(m_chance, chance_voltage_limit[i=1:numbuses], buses[i].Vmin ≤ v_response[i] ≤ buses[i].Vmax, with_probability=1-ϵ, approx="2.0")
    @constraint(m_chance, chance_θ_limit[i=1:numlines], -θᵘ ≤ θ_response[lines[i].head] - θ_response[lines[i].tail] ≤ θᵘ, with_probability=1-ϵ, approx="2.0")

    # now compute how the generators need to respond *in addition to* α*Ω
    # Δz = [ p_θv, q_pv, q_θv ]
    Δz = C*Δη

    # We must have ω_by_bus[REF] - α_by_bus[REF]Ω + ??? == Δz[1]
    # So the reference bus adjusts its production to balance the (linearized) system
    # --> ??? = Δz[1] + α_by_bus[REF]Ω - ω_by_bus[REF]
    # So the extra production at the reference bus is Δz[1] - ω_by_bus[REF]
    # (We could just set α_by_bus[REF] = 0, doesn't matter.)
    @expression(m_chance, p_response_REF, Δz[1] - sum{ω[k], k=1:numfarms; farms[k].bus==REF[1]})
    # chance constraints on active power generation
    for i in 1:numgenerators
        (generators[i].busidx != REF[1]) || continue
        @constraint(m_chance, generators[i].Pgmin <= gp_chance[i] - α_chance[i]*Ω <= generators[i].Pgmax, with_probability=1-ϵ, approx="2.0")
    end
    # if this isn't true, we have to distribute the response over the generators at the bus. Easy in principle
    @assert length(buses[REF[1]].genids) == 1
    REFGEN = buses[REF[1]].genids[1]
    @constraint(m_chance, generators[REFGEN].Pgmin ≤ gp_chance[REFGEN] + p_response_REF ≤ generators[REFGEN].Pgmax, with_probability=1-ϵ, approx="2.0")

    # q_pv and q_θv need to adjust as well
    # We must have \sum_{k at bus} γ[k]*ω[k] + ??? = Δz[...]
    for (i,busid) in enumerate(PV)
        @assert length(buses[busid].genids) == 1
        genid = buses[busid].genids[1]
        @constraint(m_chance, generators[genid].Qgmin ≤ gq_chance[genid] + Δz[1+i] - q_response[busid] ≤ generators[genid].Qgmax, with_probability=1-ϵ, approx="2.0")
    end
    gp_limits = @constraint(m_chance, generators[REFGEN].Qgmin ≤ gq_chance[REFGEN] + Δz[end] - q_response[REFGEN] ≤ generators[REFGEN].Qgmax, with_probability=1-ϵ, approx="2.0")

    # Done with the generators, now the lines

    Δx = zeros(JuMPChance.CCAffExpr, 2numbuses)
    # fill in Δη = Δ[v_pq,θ_pq,θ_pv]
    Δx[PQ] = Δη[1:length(PQ)]
    Δx[numbuses+PQ] = Δη[(length(PQ)+1):2length(PQ)]
    Δx[numbuses+PV] = Δη[(2length(PQ)+1):end]
    # Δ[pij,pji,qij,qji] ≈ vθ_to_flows_J*(Δx)
    Δflows = vθ_to_flows_J*Δx
    @assert length(Δflows) == 4numlines
    Δpij = Δflows[1:numlines]
    Δpji = Δflows[(numlines+1):2numlines]
    Δqij = Δflows[(2numlines+1):3numlines]
    Δqji = Δflows[(3numlines+1):4numlines]

    # variables to help out with the quadratic chance constraints
    @variable(m_chance, aux_pij[1:numlines] >= 0)
    @variable(m_chance, aux_pji[1:numlines] >= 0)
    @variable(m_chance, aux_qij[1:numlines] >= 0)
    @variable(m_chance, aux_qji[1:numlines] >= 0)

    @constraint(m_chance, chance_flow_limit1[i=1:numlines], -aux_pij[i] ≤ pij_chance[i] + Δpij[i] ≤ aux_pij[i], with_probability=1-ϵ, approx="1.25")
    @constraint(m_chance, chance_flow_limit2[i=1:numlines], -aux_pji[i] ≤ pji_chance[i] + Δpji[i] ≤ aux_pji[i], with_probability=1-ϵ, approx="1.25")
    @constraint(m_chance, chance_flow_limit3[i=1:numlines], -aux_qij[i] ≤ qij_chance[i] + Δqij[i] ≤ aux_qij[i], with_probability=1-ϵ, approx="1.25")
    @constraint(m_chance, chance_flow_limit4[i=1:numlines], -aux_qji[i] ≤ qji_chance[i] + Δqji[i] ≤ aux_qji[i], with_probability=1-ϵ, approx="1.25")
    @constraint(m_chance, chance_flow_limit5[i=1:numlines], aux_pij[i]^2 + aux_qij[i]^2 ≤ lines[i].u^2)
    @constraint(m_chance, chance_flow_limit6[i=1:numlines], aux_pji[i]^2 + aux_qji[i]^2 ≤ lines[i].u^2)

    status_cc = solve(m_chance, method = :Reformulate)

    println("*******************************************")
    println("Objective value: ", getobjectivevalue(m_chance))
    println("Phase angle bound (deg): ", rad2deg(θᵘ))
    println("*******************************************")


    # Define outputs of the models that are needed for calculating response functions of generators

    gp_cc   = getvalue(gp_chance)
    gq_cc   = getvalue(gq_chance)
    α_cc    = getvalue(α_chance)
    gp_det  = getvalue(gp)
    gq_det  = getvalue(gq)

    return status_det, getobjectivevalue(m), status_cc, getobjectivevalue(m_chance), getvalue(gp_chance), getvalue(gq_chance), getvalue(α_chance), getvalue(γ_chance), getvalue(r_chance), getvalue(gp), getvalue(r), getvalue(γ), total_deviation, gp_initial, gq_initial, v_initial, θ_initial,     time_det
end
