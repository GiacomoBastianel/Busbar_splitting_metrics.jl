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


gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"QCPDual" => 1, "time_limit" => 600,"MIPgap" => 5e-5, "BarHomogeneous" => 1,"BarQCPConvTol"=>1e-4)#r, "ScaleFlag"=>2, "NumericFocus"=>2,,"BarQCPConvTol"=>1e-3,) 
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 2,"linear_solver" => "ma97")
juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt, "mip_solver" => gurobi, "time_limit" => 600)
highs = JuMP.optimizer_with_attributes(HiGHS.Optimizer)#, "presolve" => 0, "presolve_dual" => 0, "simplex_dual_edge_weight_strategy" => 0, "simplex_primal_edge_weight_strategy" => 0, "simplex_price_strategy" => 0, "simplex_iteration_limit" => 100000, "simplex_perturbation" => 50.0, "simplex_dual_feasibility_tolerance" => 1.0e-9, "simplex_primal_feasibility_tolerance" => 1.0e-9, "simplex_objective_tolerance" => 1.0e-9, "simplex_infinite_bound" => 1.0e20, "simplex_infinite_cost" => 1.0e20)
##################################################################
## Processing input data
input_folder = @__DIR__

# Belgium grid without energy island
test_case_file = joinpath(dirname(dirname(input_folder)),"test_cases/pglib_opf_case39_epri.m")
test_case = _PM.parse_file(test_case_file)

s = Dict("output" => Dict("duals" => true))
result_opf = _PM.solve_opf(test_case,ACPPowerModel,ipopt,setting=s)
#result_opf_lpac = _PM.solve_opf(test_case,LPACCPowerModel,gurobi,setting=s)


powerplot(test_case)

for (b_id,b) in test_case["bus"]
    println(b_id," ",result_opf["solution"]["bus"][b_id]["lam_kcl_r"])
end

# Branches connected to 3
for (b_id,b) in test_case["branch"]
    if b["f_bus"] == 3 || b["t_bus"] == 3
        println(b_id)
    end
end

for (b_id,b) in test_case["branch"]
    println(b_id," ",result_opf["solution"]["branch"][b_id]["pf"]/b["rate_a"])
end

test_case_bs = deepcopy(test_case)

# Selecting which busbars are split -> highest dual cost
splitted_bus_ac = 3

optimizer = gurobi
formulation = LPACCPowerModel

test_case_bs,  switches_couples_ac,  extremes_ZILs_ac  = _PMTP.AC_busbar_split_AC_grid(test_case_bs,splitted_bus_ac)
test_case_bs_check = deepcopy(test_case_bs)
test_case_bs_check_auxiliary = deepcopy(test_case_bs)
results_bs = _PMTP.run_acdcsw_AC_grid_big_M(test_case_bs,formulation,optimizer)

switch_couples_feas_check = _PMTP.compute_couples_of_switches_feas_check(test_case_bs_check)
test_case_bs_check["switch_couples"] = deepcopy(switch_couples_feas_check)
test_case_bs_check_auxiliary["switch_couples"] = deepcopy(switch_couples_feas_check)

function prepare_AC_grid_feasibility_check(result_dict, input_dict, input_ac_check, switch_couples, extremes_dict,input_base)    
    orig_buses = maximum(parse.(Int, keys(input_base["bus"]))) + length(extremes_dict)
    for (sw_id,sw) in input_dict["switch"]
        if haskey(sw,"auxiliary") # Make sure ZILs are not included 
            println("sw is $(sw_id)")
            aux =  deepcopy(input_ac_check["switch"][sw_id]["auxiliary"])
            orig = deepcopy(input_ac_check["switch"][sw_id]["original"])  
            for zil in eachindex(extremes_dict)
                println("ZIL is $(zil)")
                if sw["bus_split"] == extremes_dict[zil][1] && result_dict["solution"]["switch"]["$(switch_couples[sw_id]["switch_split"])"]["status"] == 1.0  # Making sure to reconnect everything to the original if the ZIL is connected
                    if result_dict["solution"]["switch"][sw_id]["status"] >= 0.9
                        if aux == "gen"
                            input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(extremes_dict[zil][1])
                        elseif aux == "load"
                            input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(extremes_dict[zil][1])
                        elseif aux == "branch"  
                            if input_ac_check["branch"]["$(orig)"]["f_bus"] > orig_buses && switch_couples[sw_id]["bus_split"] == parse(Int64,zil)
                                input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(switch_couples[sw_id]["bus_split"])
                            elseif input_ac_check["branch"]["$(orig)"]["t_bus"] > orig_buses && switch_couples[sw_id]["bus_split"] == parse(Int64,zil)
                                input_ac_check["branch"]["$(orig)"]["t_bus"] = deepcopy(switch_couples[sw_id]["bus_split"])
                            end
                        end
                        delete!(input_ac_check["switch"],sw_id)
                    else
                        delete!(input_ac_check["switch"],sw_id)
                    end
                elseif sw["bus_split"] == extremes_dict[zil][1] && result_dict["solution"]["switch"]["$(switch_couples[sw_id]["switch_split"])"]["status"] == 0.0 # Reconnect everything to the split busbar
                    if result_dict["solution"]["switch"][sw_id]["status"] >= 0.9
                        if aux == "gen"
                            input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(sw["t_bus"])
                        elseif aux == "load"
                            input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(sw["t_bus"])
                        elseif aux == "branch" 
                            if input_ac_check["branch"]["$(orig)"]["f_bus"] > orig_buses && switch_couples[sw_id]["bus_split"] == parse(Int64,zil) 
                                    input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(sw["t_bus"])
                            elseif input_ac_check["branch"]["$(orig)"]["t_bus"] > orig_buses && switch_couples[sw_id]["bus_split"] == parse(Int64,zil)
                                if !haskey(input_ac_check["branch"]["$(orig)"],"checked")
                                    input_ac_check["branch"]["$(orig)"]["t_bus"] = deepcopy(sw["t_bus"])
                                end
                            end
                        end
                        delete!(input_ac_check["switch"],sw_id)
                    else
                        delete!(input_ac_check["switch"],sw_id)
                    end
                end
            end
        end
    end
    return input_ac_check
end
prepare_AC_grid_feasibility_check(results_bs,test_case_bs_check_auxiliary,test_case_bs_check,switch_couples_feas_check,extremes_ZILs_ac,test_case)
result_opf_check = _PMACDC.run_acdcopf(test_case_bs_check,ACPPowerModel,ipopt; setting = s)

println("Benefits from splitting busbar 3 is $((result_opf["objective"] - result_opf_check["objective"])), $((result_opf["objective"]-result_opf_check["objective"])/(result_opf["objective"])*100) %")

for (sw_id,sw) in test_case_bs_check["switch"]
    delete!(test_case_bs_check["switch"],"$sw_id")
end
function N_minus_1_checks(grid,results_dict)
    grid_base = deepcopy(grid)
    for (br_id,br) in grid["branch"]
        hourly_grid = deepcopy(grid_base)
        br["br_status"] = 0
        results_dict["$(br_id)"] = _PM.solve_opf(hourly_grid,ACPPowerModel,ipopt,setting=s)
    end
    return results_dict
end
Results_N_minus_1 = Dict{String,Any}()
@time N_minus_1_checks(test_case_bs_check,Results_N_minus_1)

powerplot(test_case_bs_check)

# Selecting which busbars are split -> highest dual cost
splitted_bus_ac_II = [3,4,14,16]

test_case_bs_II = deepcopy(test_case)

test_case_bs_II,  switches_couples_ac_II,  extremes_ZILs_ac_II  = _PMTP.AC_busbar_split_AC_grid(test_case_bs_II,splitted_bus_ac_II)
test_case_bs_check_II = deepcopy(test_case_bs_II)
test_case_bs_check_auxiliary_II = deepcopy(test_case_bs_II)
results_bs_II = _PMTP.run_acdcsw_AC_grid_big_M(test_case_bs_II,formulation,optimizer)

switch_couples_feas_check = _PMTP.compute_couples_of_switches_feas_check(test_case_bs_check_II)
test_case_bs_check_II["switch_couples"] = deepcopy(switch_couples_feas_check)
test_case_bs_check_auxiliary_II["switch_couples"] = deepcopy(switch_couples_feas_check)

prepare_AC_grid_feasibility_check(results_bs_II,test_case_bs_check_auxiliary_II,test_case_bs_check_II,switch_couples_feas_check,extremes_ZILs_ac_II,test_case)
result_opf_check_II = _PMACDC.run_acdcopf(test_case_bs_check_II,ACPPowerModel,ipopt; setting = s)

println("Benefits from splitting busbar 3 is $((result_opf["objective"] - result_opf_check_II["objective"])), $((result_opf["objective"]-result_opf_check_II["objective"])/(result_opf["objective"])*100) %")

for (sw_id,sw) in test_case_bs_check_II["switch"]
    delete!(test_case_bs_check_II["switch"],"$sw_id")
end

Results_N_minus_1_II = Dict{String,Any}()
@time N_minus_1_checks(test_case_bs_check_II,Results_N_minus_1_II)

powerplot(test_case_bs_check_II)











#Congested branches -> 3,5,7,9

##################################################################################

optimizer = gurobi
formulation = LPACCPowerModel

test_case_bs_second_stage = deepcopy(test_case_bs_check)
splitted_bus_ac_second_stage = 17
test_case_bs_second_stage,  switches_couples_ac_second_stage,  extremes_ZILs_ac_second_stage  = _PMTP.AC_busbar_split_AC_grid(test_case_bs_check,splitted_bus_ac_second_stage)
test_case_bs_check_second_stage = deepcopy(test_case_bs_check)
test_case_bs_check_auxiliary_second_stage = deepcopy(test_case_bs_check)
results_bs_second_stage = _PMTP.run_acdcsw_AC_grid_big_M(test_case_bs_second_stage,formulation,optimizer)

switch_couples_feas_check_second_stage = _PMTP.compute_couples_of_switches_feas_check(test_case_bs_check_second_stage)
test_case_bs_check_second_stage["switch_couples"] = deepcopy(switch_couples_feas_check_second_stage)
test_case_bs_check_auxiliary_second_stage["switch_couples"] = deepcopy(switch_couples_feas_check_second_stage)


prepare_AC_grid_feasibility_check(results_bs_second_stage,test_case_bs_check_auxiliary_second_stage,test_case_bs_check_second_stage,switch_couples_feas_check_second_stage,extremes_ZILs_ac_second_stage,test_case_bs_check)
result_opf_check_second_stage = _PMACDC.run_acdcopf(test_case_bs_check_second_stage,ACPPowerModel,ipopt; setting = s)
