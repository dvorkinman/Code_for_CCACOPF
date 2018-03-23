using MatpowerCases

type Bus
   nodeID::Int
   kind::Symbol
   Pd::Float64
   Qd::Float64
   gr::Float64                  # grounding resistance
   ga::Float64                  # inverse of grounding resistance
   Pg::Float64
   Qg::Float64
   Pgmax::Float64
   Qgmax::Float64
   Pgmin::Float64
   Qgmin::Float64
   pi2::Float64                 # Objective coefficient
   pi1::Float64                 # Objective coefficient
   qobjcoeff::Float64
   Pmgcost::Float64
   Vmax::Float64
   Vmin::Float64
   Ji::Float64                  # DC induced by GIC voltagex
   coord::Vector{Float64}
   genids::Vector{Int}          # Bus can contain more than one generator
   outlist::Vector{Int}         # outgoing line indices
   inlist::Vector{Int}          # incoming line indices
   Gs::Float64     #shunt conductance
   Bs::Float64     #shunt susceptance
   Vm::Float64     #voltage magnitude (for PV buses)
   function Bus(nodeID, kind, Pd, Qd, Vmax, Vmin, Gs, Bs, Vm)
      if kind == 1
          k = :PQ
      elseif kind == 2
          k = :PV
      elseif kind == 3
          k = :Ref
      else
          error()
      end
      b = new(nodeID, k, Pd, Qd)
      b.gr = 0
      b.ga = 0
      b.Pg = 0
      b.Qg = 0
      b.Pgmax = 0
      b.Qgmax = 0
      b.pi1 = 0
      b.pi2 = 0
      b.Pmgcost = 0.0
      b.Vmax = Vmax
      b.Vmin = Vmin
      b.Ji = 0.0
      b.coord=[0.0, 0.0]
      b.genids = Int[]
      b.outlist = Int[]
      b.inlist = Int[]
      b.Gs = Gs
      b.Bs = Bs
      b.Vm = Vm
      return b
   end
end

# Functions for bus
function setg(b::Bus, genidx, Pg, Qg, Pgmax, Pgmin, Qgmax, Qgmin)
   b.Pg += Pg
   b.Qg += Qg
   b.Pgmax += Pgmax
   b.Pgmin += Pgmin
   b.Qgmax += Qgmax
   b.Qgmin += Qgmin
   if b.kind == :PQ
      warn("Generator $genidx was assigned to bus $(b.nodeID), but this bus has type PV")
   end
   push!(b.genids,genidx)
end

type Generator
   genID::Int
   busidx::Int
   Pg::Float64
   Qg::Float64
   Pgmax::Float64
   Pgmin::Float64
   Qgmax::Float64
   Qgmin::Float64
   pi1::Float64
   pi2::Float64
   pi3::Float64
   function Generator(genID, busidx, Pg, Qg, Pgmax, Pgmin, Qgmax, Qgmin)
      g = new(genID, busidx, Pg, Qg, Pgmax, Pgmin, Qgmax, Qgmin)
      g.pi1 = 1.0
      g.pi2 = 1.0
      g.pi3 = 1.0
      return g
   end
end

type Line
   arcID::Int
   tail::Int # the "to" node
   head::Int # the "from" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   γ::Float64 # Conductance
   β::Float64 # Susceptance
   u::Float64 # the capacity of the line
   ratio::Float64 # turns ratio of transformer on the line
   distance_scale::Float64 # this will be used to scale u
   Imax::Float64 # max Iᵃ
   Imin::Float64 # min Iᵃ
   a::Float64 # inverse of line resistance
   Jij::Float64 # induced DC flow by GMDs
   Iij::Float64 # current causing magnetic saturation of transformer
   b_charge:: Float64 # charging susceptance
   function Line(arcID, tail, head, r, x, u, turns, d, Imax, Imin,b_charge)
      line = new(arcID, tail, head, r, x)
      line.γ = r/(r^2+x^2)
      line.β = -x/(r^2+x^2)
      line.u = u
      line.ratio = turns
      line.distance_scale = d
      line.Imax = Imax
      line.Imin = Imin
      line.a= 1/r
      line.Jij= 0
      line.Iij= Imax
      line.b_charge = b_charge
      return line
   end
end

getThermalCapacity(l::Line, mvaBase) = l.u#/mvaBase  # line limits
getSyncCapacity(l::Line, mvaBase) = l.y

type TransLine
   translineID::Int
   arcID::Int
   tail::Int # the "to" node
   head::Int # the "from" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   γ::Float64 # Conductance
   β::Float64 # Susceptance
   u::Float64 # the capacity of the line
   ratio::Float64 # turns ratio of the line
   distance_scale::Float64 # this will be used to scale u
   Imax::Float64 # max Iᵃ
   Imin::Float64 # min Iᵃ
   a::Float64 # inverse of line resistance
   Jij::Float64 # induced DC flow by GMDs
   Iij::Float64 # current causing magnetic saturation of transformer
   b_charge::Float64
   function TransLine(translineID, arcID, tail, head, r, x, u, turns, d, Imax, Imin, b_charge)
      transline = new(translineID, arcID, tail, head, r, x)
      transline.γ = r/(r^2+x^2)
      transline.β = -x/(r^2+x^2)
      transline.u = u
      transline.ratio = turns
      transline.distance_scale = d
      transline.Imax = Imax
      transline.Imin = Imin
      transline.a= 1/r
      transline.Jij= 0.0
      transline.Iij= Imax
      transline.b_charge=b_charge
      return transline
   end
end

type Scenario
   lineIDs::Vector{Int}
   function Scenario(lineIDs)
      s = new(lineIDs)
      return s
   end
end

function GetlineID(lines, head, tail)
   lineIDs = Int[]
   for i=1:length(lines)
      if((lines[i].head == head && lines[i].tail == tail) || (lines[i].head == tail && lines[i].tail == head))
         push!(lineIDs,lines[i].arcID)  #Push all the repeated lines between given nodes.
      end
   end
   return lineIDs
end

type Farm
    μ::Float64
    σ::Float64
    bus::Int
end

function readcase(casefilename, extras, loadscale, thermalLimitscale, mvaBase)

   case = loadcase(convert(ASCIIString, casefilename))

   ## bus data
   # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
   busmat = case["bus"]
   busIDmap = Dict()

   buses = Bus[]
   for i in 1:size(busmat,1)
      nodeID = i
      busIDmap[busmat[i,1]] = i
      bustype = round(Int, busmat[i,2])
      Pd = busmat[i,3]*loadscale
      Qd = busmat[i,4]*loadscale
      Gs = busmat[i,5]
      if busmat[i,6] < 0
        Bs= busmat[i,6]
      else
        Bs = busmat[i,6]
      end
      area = busmat[i,7]
      Vm = busmat[i,8]
      Va = busmat[i,9]
      baseKV = busmat[i,10]
      zone = busmat[i,11]
      Vmax = busmat[i,12]
      Vmin = busmat[i,13]
      b = Bus(nodeID, bustype, Pd./mvaBase, Qd./mvaBase, Vmax, Vmin, Gs/mvaBase, Bs/mvaBase, Vm)
      push!(buses, b)
   end
   ####################################################################################################################

   ## generator data
   # bus Pg Qg	Qmax Qmin Vg mBase status Pmax Pmin Pc1 Pc2	Qc1min Qc1max Qc2min Qc2max ramp_agc ramp_10 ramp_30 ramp_q apf
   generatorlist = Int[]
   generators = Generator[]
   genmat = case["gen"]
   gencounter = 1
   for i in 1:size(genmat,1)
      busidx = busIDmap[genmat[i,1]]    #indicate generator i is connected to bus = busidx
      Pg = genmat[i,2]
      Qg = genmat[i,3]
      Qgmax = genmat[i,4]
      Qgmin = genmat[i,5]
      Pgmax = genmat[i,9]
      Pgmin = genmat[i,10]
      g = Generator(gencounter, busidx, Pg/mvaBase, Qg/mvaBase, Pgmax/mvaBase, Pgmin/mvaBase, Qgmax/mvaBase, Qgmin/mvaBase)
      push!(generators, g)
      push!(generatorlist, busidx)
      setg(buses[busidx], gencounter, g.Pg, g.Qg, g.Pgmax, g.Pgmin, g.Qgmax, g.Qgmin)
      gencounter = gencounter + 1
   end

   for i in 1:length(buses)
      if buses[i].kind == 2 && length(buses[i].genids) == 0
         warn("Bus $i is listed as a generator, but no corresponding generator information found, changing to kind 1")
         buses[i].kind = 1
      end
   end

   gencost = case["gencost"]
   for g in 1:length(generators)
      generators[g].pi1 = gencost[g,5]*mvaBase^2
      generators[g].pi2 = gencost[g,6]*mvaBase
      generators[g].pi3 = gencost[g,7]
   end

   ######################################################################################################
   ## branch data
   # fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax
   branchmat = case["branch"]
   lines = Line[]
   translines = TransLine[]
   translineID = 0
   for i in 1:size(branchmat,1)
      fbus = busIDmap[branchmat[i,1]]
      tbus = busIDmap[branchmat[i,2]]
      x = branchmat[i,4]
      r = branchmat[i,3]
      u = branchmat[i,6]/mvaBase
      if branchmat[i,9]!=0
        turns = branchmat[i,9]
      else
        turns = 1.0
      end
      b_charge = branchmat[i,5]
      Imax = u*thermalLimitscale/buses[fbus].Vmin
      Imin = 0
      push!(buses[fbus].outlist, i)
      push!(buses[tbus].inlist, i)
      l = Line(i, tbus, fbus, r, x, u*thermalLimitscale, turns, thermalLimitscale, Imax, Imin, b_charge)  # turns==1 if no transformer on the line
      push!(lines,l)
      # if this line has a transformer, then add this line to subset translines
      # and set up the grouding resistance of head and tail bus of the line.
      if turns != 1
         translineID += 1
#         println("Line $i has a transformer")
#         println("transformer is sitting at $(lines[i].head)")
         tl = TransLine(translineID, i, tbus, fbus, r, x, u*thermalLimitscale, turns, thermalLimitscale, Imax, Imin, b_charge) # i is translineID
         push!(translines, tl)
         buses[fbus].ga = 0.1  # transformers are siting at the heads of lines
         #buses[fbus].ga = setga(buses[fbus], r)   # transformers are siting at the heads of lines
     end
   end

   ### Compute Jij
   for i=1:length(lines)
      from=lines[i].head
      to= lines[i].tail
      lines[i].Jij=(buses[from].Ji-buses[to].Ji) # J_ij= J_i - J_j
   end

   generatorlist = unique(generatorlist)
   numbuses = length(buses)
   return generators, generatorlist, buses, lines, translines
end

function readconfig(configfilename)
   println("\nreading config $configfilename")
   refbus = 0
   lines = readlines(open(configfilename,"r"))
   numlines = length(lines)

   logfilename = "none"
   casefilename = "none"
   lines_cost_filename = "none"

   lookingforend = true
   mvaBase = 100
   loadscale = 1.2
   thermalLimitscale = 0.8
   extras = Dict()

   for l in lines
      startswith(l,'#') && continue

      thisline = split(l)
      length(thisline) > 0 || continue
      if thisline[1] == "END"
         break
      elseif thisline[1] == "case"
         casefilename = thisline[2]
      elseif thisline[1] == "refbus"
         refbus = parse(Int,thisline[2])
      elseif thisline[1] == "mvaBase"
         mvaBase = float(thisline[2])
         println(">>>> mvaBase = $mvaBase")
      elseif thisline[1] == "loadscale"
         loadscale = float(thisline[2])
         println(">>>> loadscale = $loadscale")
      elseif thisline[1] == "thermalLimitscale"
         thermalLimitscale = float(thisline[2])
         println(">>>> thermalLimitscale = $thermalLimitscale")
      elseif thisline[1] == "logfile"
         logfilename = thisline[2]
         println(">>>> logfilename = $logfilename")
      else
         extras[thisline[1]] = thisline[2]
         if thisline[1] == "theta_u"
            println(">>>> $(thisline[1]) = $(thisline[2])")
         elseif thisline[1] =="lambda"
            println(">>>> $(thisline[1]) = $(thisline[2])")
         end
      end
   end

   return casefilename, refbus, loadscale, thermalLimitscale, extras, mvaBase
end
