using PowerModels; const _PM = PowerModels
using PowerModelsACDC; const _PMACDC = PowerModelsACDC
using PowerModelsTopologicalActionsII; const _PMTP = PowerModelsTopologicalActionsII
using DIRECTIONS_WP4; const _WP4 = DIRECTIONS_WP4 
using Gurobi
using JuMP
using HiGHS
using DataFrames
using CSV
using Feather
using JSON
using Ipopt
using Plots
using StatsBase
using Juniper
using Statistics
using HSL_jll
using PowerPlots


gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"QCPDual" => 1, "time_limit" => 600,"MIPgap" => 2.5e-3, "BarHomogeneous" => 1,"BarQCPConvTol"=>1e-4)#r, "ScaleFlag"=>2, "NumericFocus"=>2,,"BarQCPConvTol"=>1e-3,) 
gurobi_no_dual = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "time_limit" => 600,"MIPgap" => 1e-2, "BarHomogeneous" => 1)#r, "ScaleFlag"=>2, "NumericFocus"=>2,,"BarQCPConvTol"=>1e-3,) 

ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 2,"linear_solver" => "ma97")
juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt, "mip_solver" => gurobi, "time_limit" => 600)
highs = JuMP.optimizer_with_attributes(HiGHS.Optimizer)#, "presolve" => 0, "presolve_dual" => 0, "simplex_dual_edge_weight_strategy" => 0, "simplex_primal_edge_weight_strategy" => 0, "simplex_price_strategy" => 0, "simplex_iteration_limit" => 100000, "simplex_perturbation" => 50.0, "simplex_dual_feasibility_tolerance" => 1.0e-9, "simplex_primal_feasibility_tolerance" => 1.0e-9, "simplex_objective_tolerance" => 1.0e-9, "simplex_infinite_bound" => 1.0e20, "simplex_infinite_cost" => 1.0e20)
##################################################################
## Processing input data
input_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/data_sources"

# Belgium grid without energy island
test_case_file = joinpath(input_folder,"CaliforniaTestSystem.m")
test_case = _PM.parse_file(test_case_file)

s = Dict("output" => Dict("duals" => true))
s_no_duals = Dict("output" => Dict("duals" => false))

result_opf = _PM.solve_opf(test_case,ACPPowerModel,ipopt,setting=s)
result_opf_lpac = _PM.solve_opf(test_case,LPACCPowerModel,gurobi_no_dual,setting=s_no_duals)


powerplot(test_case)

for (b_id,b) in test_case["bus"]
    println(b_id," ",result_opf["solution"]["bus"][b_id]["lam_kcl_r"])
end
duals = []
for (b_id,b) in test_case["bus"]
    push!(duals,[parse(Int64,b_id),abs(result_opf["solution"]["bus"][b_id]["lam_kcl_r"])])
end

sorted_duals = sort(duals, by = x -> x[2], rev = true)


for (b_id,b) in test_case["branch"]
    println(b_id," ",abs(result_opf["solution"]["branch"][b_id]["pf"]/b["rate_a"]))
end

for (b_id,b) in test_case["branch"]
    println(b_id," f_bus ",b["f_bus"]," t_bus ",b["t_bus"]," ",abs(abs(result_opf["solution"]["bus"]["$(b["f_bus"])"]["lam_kcl_r"]) - abs(result_opf["solution"]["bus"]["$(b["t_bus"])"]["lam_kcl_r"])))
end
duals_diff = []
for (b_id,b) in test_case["branch"]
    push!(duals_diff,[parse(Int64,b_id),abs(abs(result_opf["solution"]["bus"]["$(b["f_bus"])"]["lam_kcl_r"]) - abs(result_opf["solution"]["bus"]["$(b["t_bus"])"]["lam_kcl_r"]))])
end
sorted_duals_diff = sort(duals_diff, by = x -> x[2], rev = true)

buses_split = []
for i in 1:50
    push!(buses_split,Int64(duals_diff[i][1]))
end

test_case_bs_try = deepcopy(test_case)
optimizer = gurobi
formulation = LPACCPowerModel

result = Dict{String,Any}()
result_bs = Dict{String,Any}()
for b_id in buses_split
    splitted_bus_ac = b_id
    println("----------------")
    println("Bus")
    println(splitted_bus_ac)
    println("----------------")
    hourly_test_case = deepcopy(test_case_bs_try)
    hourly_test_case,  switches_couples_ac,  extremes_ZILs_ac  = _PMTP.AC_busbar_split_AC_grid(hourly_test_case,splitted_bus_ac)
    hourly_test_case_bs_check = deepcopy(hourly_test_case)
    hourly_test_case_bs_check_auxiliary = deepcopy(hourly_test_case)
    println("----------------")
    println("Switch couples")
    println(switches_couples_ac)
    println("----------------")
    results_bs_hourly = _PMTP.run_acdcsw_AC_big_M_ZIL(hourly_test_case,formulation,gurobi_no_dual)
    result_bs["$b_id"] = deepcopy(results_bs_hourly)
    #switch_couples_feas_check = _PMTP.compute_couples_of_switches_feas_check(hourly_test_case_bs_check)
    #hourly_test_case_bs_check["switch_couples"] = deepcopy(switch_couples_feas_check)
    #hourly_test_case_bs_check_auxiliary["switch_couples"] = deepcopy(switch_couples_feas_check)
    #prepare_AC_grid_feasibility_check(results_bs_hourly,hourly_test_case_bs_check_auxiliary,hourly_test_case_bs_check,switch_couples_feas_check,extremes_ZILs_ac,hourly_test_case)
    #result_opf_check = _PMACDC.run_acdcopf(hourly_test_case_bs_check,ACPPowerModel,ipopt; setting = s)
    #result["$b_id"] = deepcopy(result_opf_check)
end

bs_buses = []
println("Start optimization")
for i in keys(result)
    if result[i]["objective"] < result_opf["objective"]
        push!(bs_buses,i)        
        println("Busbar $(i) $(result[i]["objective"])")
        println("----------------")
        println("Reduction in generation cost $((1 - result[i]["objective"]/result_opf["objective"])*100)%")
        println(" ")
    end
end

for (br_id,br) in test_case_bs_try["branch"]
    if br["f_bus"] == 308 || br["t_bus"] == 308
        println("Branch $(br_id) $(br["f_bus"]) $(br["t_bus"])")
    end
end

diff_opf_bs = []




obj_bs = [result_bs[b_id]["objective"] for (b_id,b) in test_case_bs_try["bus"]]

plot(obj_bs)

buses = []
bs_opf_diff = []
dual_buses = []
for (b_id,b) in test_case_bs_try["bus"]
    push!(buses,parse(Int64,b_id))
    push!(dual_buses,abs(result_opf["solution"]["bus"][b_id]["lam_kcl_r"]))
    push!(bs_opf_diff,result_opf["objective"] - result[b_id]["objective"])
end
zero = 0.0*collect(1:length(bs_opf_diff))
vertical_line = 100*collect(1:length(bs_opf_diff))
scatter(bs_opf_diff/10^3,dual_buses,legend= :none,grid = :none,ylabel = "Abs duals lam_kcl_r, AC-feasibility check",xlabel = "Diff AC-OPF - AC-feasibility check [k€], > 0 = reduction in costs",ylims = (3000,5000),xlims = (-30,15),title = "RTS GMC test case, busbar splitting, no OTS")
annotate!([(bs_opf_diff[i]/10^3, dual_buses[i], text(buses[i], 9, :red, :right)) for i in 1:length(buses) if bs_opf_diff[i] > 0])
plot!(zero,vertical_line,legend= :none,color = :red)





