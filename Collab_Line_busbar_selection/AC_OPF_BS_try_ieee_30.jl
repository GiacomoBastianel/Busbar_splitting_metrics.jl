using PowerModels; const _PM = PowerModels
using PowerModelsACDC; const _PMACDC = PowerModelsACDC
using PowerModelsTopologicalActionsII; const _PMTP = PowerModelsTopologicalActionsII
using Gurobi, JuMP, HiGHS, DataFrames, CSV, Feather, JSON, Ipopt, Plots, StatsBase, Juniper, Statistics, PowerPlots


gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "time_limit" => 600,"MIPgap" => 5e-5, "BarHomogeneous" => 1)#,"QCPDual" => 1,"BarQCPConvTol"=>1e-4)#r, "ScaleFlag"=>2, "NumericFocus"=>2,,"BarQCPConvTol"=>1e-3,) 
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 2)#,"linear_solver" => "ma97")
juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt, "mip_solver" => gurobi, "time_limit" => 600)
highs = JuMP.optimizer_with_attributes(HiGHS.Optimizer)#, "presolve" => 0, "presolve_dual" => 0, "simplex_dual_edge_weight_strategy" => 0, "simplex_primal_edge_weight_strategy" => 0, "simplex_price_strategy" => 0, "simplex_iteration_limit" => 100000, "simplex_perturbation" => 50.0, "simplex_dual_feasibility_tolerance" => 1.0e-9, "simplex_primal_feasibility_tolerance" => 1.0e-9, "simplex_objective_tolerance" => 1.0e-9, "simplex_infinite_bound" => 1.0e20, "simplex_infinite_cost" => 1.0e20)
##################################################################
## Processing input data
input_folder = @__DIR__

# Belgium grid without energy island
test_case_file = joinpath(dirname(input_folder),"test_cases/pglib_opf_case30_ieee.m")
test_case = _PM.parse_file(test_case_file)

s_dual = Dict("output" => Dict("duals" => true))
s = Dict("output" => Dict("duals" => false))
result_opf = _PM.solve_opf(test_case,ACPPowerModel,ipopt,setting=s_dual)
result_opf_lpac = _PM.solve_opf(test_case,LPACCPowerModel,gurobi,setting=s)


for (b_id,b) in test_case["bus"]
    println(b_id," ",result_opf["solution"]["bus"][b_id]["lam_kcl_r"])
end
duals = []
for (b_id,b) in test_case["bus"]
    push!(duals,[parse(Int64,b_id),result_opf["solution"]["bus"][b_id]["lam_kcl_r"]])
end
sorted_duals = sort(duals, by = x -> x[2], rev = true)

test_case_bs = deepcopy(test_case)
test_case_bs_example = deepcopy(test_case)

###############################################################################
# Splitting busbars
function split_one_bus_per_time(test_case,results_dict,results_dict_ac_check,results_dict_lpac_check)
    for (b_id,b) in test_case["bus"]
        results_dict["$b_id"] = Dict{String,Any}()
        results_dict_ac_check["$b_id"] = Dict{String,Any}()
        test_case_bs = deepcopy(test_case)
        splitted_bus_ac = parse(Int64,b_id)
        test_case_bs,  switches_couples_ac,  extremes_ZILs_ac  = _PMTP.AC_busbar_split_AC_grid(test_case_bs,splitted_bus_ac)
        test_case_bs["switch"]["1"]["cost"] = 50.0
        results_dict["$b_id"] = _PMTP.run_acdcsw_AC_grid_big_M(test_case_bs,LPACCPowerModel,gurobi)
        test_case_bs_check = deepcopy(test_case_bs)
        test_case_bs_check_auxiliary = deepcopy(test_case_bs)
        if results_dict["$b_id"]["termination_status"] == JuMP.OPTIMAL
            prepare_AC_grid_feasibility_check(results_dict["$b_id"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples_ac,extremes_ZILs_ac,test_case)
            results_dict_ac_check["$b_id"] = deepcopy(_PMACDC.run_acdcopf(test_case_bs_check,ACPPowerModel,ipopt; setting = s))
            results_dict_lpac_check["$b_id"] = deepcopy(_PMACDC.run_acdcopf(test_case_bs_check,LPACCPowerModel,gurobi; setting = s))
        end
    end
end

function busbar_splitting_tree(test_case,results_dict,results_dict_ac_check,opf_lpac)
    for (b_id,b) in test_case["bus"]
        results_dict["$b_id"] = Dict{String,Any}()
        results_dict_ac_check["$b_id"] = Dict{String,Any}()
        test_case_bs = deepcopy(test_case)
        splitted_bus_ac = parse(Int64,b_id)
        test_case_bs,  switches_couples_ac,  extremes_ZILs_ac  = _PMTP.AC_busbar_split_AC_grid(test_case_bs,splitted_bus_ac)
        results_dict["$b_id"] = _PMTP.run_acdcsw_AC_grid_big_M(test_case_bs,LPACCPowerModel,gurobi)
        test_case_bs_check = deepcopy(test_case_bs)
        test_case_bs_check_auxiliary = deepcopy(test_case_bs)
        if results_dict["$b_id"]["termination_status"] == JuMP.OPTIMAL
            prepare_AC_grid_feasibility_check(results_dict["$b_id"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples_ac,extremes_ZILs_ac,test_case)
            results_dict_ac_check["$b_id"] = _PMACDC.run_acdcopf(test_case_bs_check,LPACCPowerModel,ipopt; setting = s)
            if results_dict["$b_id"]["objective"] < opf_lpac
                for (b_2_id,b_2) in test_case["bus"]
                    if b_2_id != b_id
                        test_case_bs = deepcopy(test_case)
                        splitted_bus_ac = [b_id, parse(Int64,b_2_id)]
                        test_case_bs,  switches_couples_ac,  extremes_ZILs_ac  = _PMTP.AC_busbar_split_AC_grid(test_case_bs,splitted_bus_ac)
                        results_b_2 = _PMTP.run_acdcsw_AC_grid_big_M(test_case_bs,LPACCPowerModel,gurobi)
                        if results_b_2["termination_status"] == JuMP.OPTIMAL && results_b_2["objective"] < results_dict["$b_id"]["objective"]
                            results_dict["$b_id"]["$b_2_id"] = Dict{String,Any}()
                            results_dict_ac_check["$b_id"]["$b_2_id"] = Dict{String,Any}()
                            results_dict["$b_id"]["$b_2_id"] = deepcopy(results_b_2)
                            test_case_bs_check = deepcopy(test_case_bs)
                            test_case_bs_check_auxiliary = deepcopy(test_case_bs)
                            prepare_AC_grid_feasibility_check(results_dict["$b_id"]["$b_2_id"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples_ac,extremes_ZILs_ac,test_case)
                            results_dict_ac_check["$b_id"]["$b_2_id"] = _PMACDC.run_acdcopf(test_case_bs_check,LPACCPowerModel,ipopt; setting = s)
                        end
                    end
                end
            end
        end
    end
end

function N_minus_1_checks(grid,results_dict,formulation,optimizer)
    grid_base = deepcopy(grid)
    for (br_id,br) in grid["branch"]
        hourly_grid = deepcopy(grid_base)
        br["br_status"] = 0
        hourly_result = _PM.solve_opf(hourly_grid,formulation,optimizer,setting=s)
        if hourly_result["termination_status"] != JuMP.OPTIMAL
            results_dict["$(br_id)"] = deepcopy(hourly_result)
        end
    end
    return results_dict
end

function split_one_bus_per_time_with_N_1_check(test_case,results_dict,results_dict_ac_check,results_dict_lpac_check)
    for (b_id,b) in test_case["bus"]
        results_dict["$b_id"] = Dict{String,Any}()
        results_dict_ac_check["$b_id"] = Dict{String,Any}()
        test_case_bs = deepcopy(test_case)
        splitted_bus_ac = parse(Int64,b_id)
        test_case_bs,  switches_couples_ac,  extremes_ZILs_ac  = _PMTP.AC_busbar_split_AC_grid(test_case_bs,splitted_bus_ac)
        test_case_bs["switch"]["1"]["cost"] = 50.0
        results_dict["$b_id"] = _PMTP.run_acdcsw_AC_grid_big_M(test_case_bs,LPACCPowerModel,gurobi)
        test_case_bs_check = deepcopy(test_case_bs)
        test_case_bs_check_auxiliary = deepcopy(test_case_bs)
        if results_dict["$b_id"]["termination_status"] == JuMP.OPTIMAL
            prepare_AC_grid_feasibility_check(results_dict["$b_id"],test_case_bs_check_auxiliary,test_case_bs_check,switches_couples_ac,extremes_ZILs_ac,test_case)
            results_dict_ac_check["$b_id"] = _PMACDC.run_acdcopf(test_case_bs_check,ACPPowerModel,ipopt; setting = s)
            results_dict_lpac_check["$b_id"] = _PMACDC.run_acdcopf(test_case_bs_check,LPACCPowerModel,gurobi; setting = s)
            results_dict_lpac_check["$b_id"]["feasibility_check"] = Dict{String,Any}()
            results_dict_ac_check["$b_id"]["feasibility_check"] = Dict{String,Any}()
            N_minus_1_checks(test_case_bs_check,results_dict_lpac_check["$b_id"]["feasibility_check"],LPACCPowerModel,gurobi)
            N_minus_1_checks(test_case_bs_check,results_dict_ac_check["$b_id"]["feasibility_check"],ACPPowerModel,ipopt)
        end
    end
end

function prepare_AC_grid_feasibility_check(result_dict, input_dict, input_ac_check, switch_couples, extremes_dict, input_base)
    orig_buses = maximum(parse.(Int, keys(input_base["bus"])))
    println("Number of original buses is $orig_buses")
    for (sw_id,sw) in input_dict["switch"]
            if !haskey(sw,"auxiliary")
                println("SWITCH $sw_id, BUS $(sw["t_bus"])")
                if result_dict["solution"]["switch"][sw_id]["status"] >= 0.9
                    println("Switch $sw_id is closed, Connecting everything back, no busbar splitting on bus $(sw["bus_split"])")
                    #delete!(input_ac_check["bus"],"$(input_ac_check["switch"][sw_id]["t_bus"])")
                    for l in keys(switch_couples)
                        if switch_couples[l]["bus_split"] == sw["bus_split"]
                            println("SWITCH COUPLE IS $l")
                            #println("Starting from switch $(switch_couples[l]["f_sw"]), with t_bus $(input_ac_check["switch"]["$(switch_couples[l]["f_sw"])"]["t_bus"])")
                            #println("Then switch $(switch_couples[l]["t_sw"]), with t_bus $(input_ac_check["switch"]["$(switch_couples[l]["t_sw"])"]["t_bus"])")
                            if input_dict["switch"]["$(switch_couples[l]["f_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                println("WE GO WITH SWITCH $(switch_couples[l]["f_sw"])")
                                aux =  deepcopy(input_dict["switch"]["$(switch_couples[l]["f_sw"])"]["auxiliary"])
                                orig = deepcopy(input_dict["switch"]["$(switch_couples[l]["f_sw"])"]["original"])
                                if aux == "gen"
                                    delete!(input_ac_check["bus"],input_ac_check["gen"]["$(orig)"]["gen_bus"])
                                    input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif aux == "load"
                                    delete!(input_ac_check["bus"],input_ac_check["load"]["$(orig)"]["load_bus"])
                                    input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif aux == "convdc"
                                    delete!(input_ac_check["bus"],input_ac_check["convdc"]["$(orig)"]["busac_i"])
                                    input_ac_check["convdc"]["$(orig)"]["busac_i"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif aux == "branch"                
                                    if input_dict["branch"]["$(orig)"]["f_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["f_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                        delete!(input_ac_check["bus"],input_ac_check["branch"]["$(orig)"]["f_bus"])
                                        input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                    elseif input_dict["branch"]["$(orig)"]["t_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["f_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                        delete!(input_ac_check["bus"],input_ac_check["branch"]["$(orig)"]["t_bus"])
                                        input_ac_check["branch"]["$(orig)"]["t_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                    end
                                end
                            elseif input_dict["switch"]["$(switch_couples[l]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                println("WE GO WITH SWITCH $(switch_couples[l]["t_sw"])")
                                println("----------------------------")
                                aux =  deepcopy(input_dict["switch"]["$(switch_couples[l]["t_sw"])"]["auxiliary"])
                                orig = deepcopy(input_dict["switch"]["$(switch_couples[l]["t_sw"])"]["original"])
                                if aux == "gen"
                                    delete!(input_ac_check["bus"],input_ac_check["gen"]["$(orig)"]["gen_bus"])
                                    input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif aux == "load"
                                    delete!(input_ac_check["bus"],input_ac_check["load"]["$(orig)"]["load_bus"])
                                    input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif aux == "convdc"
                                    delete!(input_ac_check["bus"],input_ac_check["convdc"]["$(orig)"]["busac_i"])
                                    input_ac_check["convdc"]["$(orig)"]["busac_i"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif aux == "branch"                
                                    if input_dict["branch"]["$(orig)"]["f_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                        delete!(input_ac_check["bus"],input_ac_check["branch"]["$(orig)"]["f_bus"])
                                        input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                    elseif input_dict["branch"]["$(orig)"]["t_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                        delete!(input_ac_check["bus"],input_ac_check["branch"]["$(orig)"]["t_bus"])
                                        input_ac_check["branch"]["$(orig)"]["t_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                    end
                                end
                            end
                        end
                    end
                elseif result_dict["solution"]["switch"][sw_id]["status"] <= 0.1
                    println("Switch $sw_id is open, busbar splitting on bus $(sw["bus_split"])")
                    delete!(input_ac_check["switch"],sw_id)
                    for l in keys(switch_couples)
                        if switch_couples[l]["bus_split"] == sw["bus_split"]
                            println("SWITCH COUPLE IS $l")
                            #println("Starting from switch $(switch_couples[l]["f_sw"]), with t_bus $(input_ac_check["switch"]["$(switch_couples[l]["f_sw"])"]["t_bus"])")
                            #println("Then switch $(switch_couples[l]["t_sw"]), with t_bus $(input_ac_check["switch"]["$(switch_couples[l]["t_sw"])"]["t_bus"])")
                            if result_dict["solution"]["switch"]["$(switch_couples["$l"]["switch_split"])"]["status"] >= 0.9
                                aux =  deepcopy(input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["auxiliary"])
                                orig = deepcopy(input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["original"])
                                if aux == "gen"
                                    delete!(input_ac_check["bus"],input_ac_check["gen"]["$(orig)"]["gen_bus"])
                                    input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif aux == "load"
                                    delete!(input_ac_check["bus"],input_ac_check["load"]["$(orig)"]["load_bus"])
                                    input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif aux == "convdc"
                                    delete!(input_ac_check["bus"],input_ac_check["convdc"]["$(orig)"]["busac_i"])
                                    input_ac_check["convdc"]["$(orig)"]["busac_i"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif aux == "branch"                
                                    if input_dict["branch"]["$(orig)"]["f_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                        delete!(input_ac_check["bus"],input_ac_check["branch"]["$(orig)"]["f_bus"])
                                        input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                    elseif input_dict["branch"]["$(orig)"]["t_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                        delete!(input_ac_check["bus"],input_ac_check["branch"]["$(orig)"]["t_bus"])
                                        input_ac_check["branch"]["$(orig)"]["t_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                    end
                                end
                            elseif result_dict["solution"]["switch"]["$(switch_couples["$l"]["switch_split"])"]["status"] <= 0.1
                                switch_t = deepcopy(input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]) 
                                aux_t = switch_t["auxiliary"]
                                orig_t = switch_t["original"]
                                print([l,aux_t,orig_t],"\n")
                                if result_dict["solution"]["switch"]["$(switch_t["index"])"]["status"] <= 0.1
                                    delete!(input_ac_check["switch"],"$(switch_t["index"])")
                                elseif result_dict["solution"]["switch"]["$(switch_t["index"])"]["status"] >= 0.9
                                    if aux_t == "gen"
                                        delete!(input_ac_check["bus"],input_ac_check["gen"]["$(orig_t)"]["gen_bus"])
                                        input_ac_check["gen"]["$(orig_t)"]["gen_bus"] = deepcopy(input_dict["switch"]["$(switch_t["index"])"]["t_bus"]) # here it needs to be the bus of the switch
                                        print([l,aux_t,orig_t,input_ac_check["gen"]["$(orig_t)"]["gen_bus"]],"\n")
                                    elseif aux_t == "load"
                                        delete!(input_ac_check["bus"],input_ac_check["load"]["$(orig_t)"]["load_bus"])
                                        input_ac_check["load"]["$(orig_t)"]["load_bus"] = deepcopy(input_dict["switch"]["$(switch_t["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["load"]["$(orig_t)"]["load_bus"]],"\n")
                                    elseif aux_t == "convdc"
                                        delete!(input_ac_check["bus"],input_ac_check["convdc"]["$(orig_t)"]["busac_i"])
                                        input_ac_check["convdc"]["$(orig_t)"]["busac_i"] = deepcopy(input_dict["switch"]["$(switch_t["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["convdc"]["$(orig_t)"]["busac_i"]],"\n")
                                    elseif aux_t == "branch" 
                                        if input_ac_check["branch"]["$(orig_t)"]["f_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
                                            delete!(input_ac_check["bus"],input_ac_check["branch"]["$(orig_t)"]["f_bus"])
                                            input_ac_check["branch"]["$(orig_t)"]["f_bus"] = deepcopy(input_dict["switch"]["$(switch_t["index"])"]["t_bus"])
                                            print([l,aux_t,orig_t,input_ac_check["branch"]["$(orig_t)"]["f_bus"]],"\n")
                                        elseif input_ac_check["branch"]["$(orig_t)"]["t_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
                                            delete!(input_ac_check["bus"],input_ac_check["branch"]["$(orig_t)"]["t_bus"])
                                            input_ac_check["branch"]["$(orig_t)"]["t_bus"] = deepcopy(input_dict["switch"]["$(switch_t["index"])"]["t_bus"])
                                            print([l,aux_t,orig_t,input_ac_check["branch"]["$(orig_t)"]["t_bus"]],"\n")
                                        end
                                    end
                                end
                            
                                switch_f = deepcopy(input_dict["switch"]["$(switch_couples["$l"]["f_sw"])"]) 
                                aux_f = switch_f["auxiliary"]
                                orig_f = switch_f["original"]
                                print([l,aux_f,orig_f],"\n")
                                if result_dict["solution"]["switch"]["$(switch_f["index"])"]["status"] <= 0.1
                                    delete!(input_ac_check["switch"],"$(switch_t["index"])")
                                else
                                    if aux_f == "gen"
                                        delete!(input_ac_check["bus"],input_ac_check["gen"]["$(orig_f)"]["gen_bus"])
                                        input_ac_check["gen"]["$(orig_f)"]["gen_bus"] = deepcopy(input_dict["switch"]["$(switch_f["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["gen"]["$(orig_f)"]["gen_bus"]],"\n")
                                    elseif aux_f == "load"
                                        delete!(input_ac_check["bus"],input_ac_check["load"]["$(orig_f)"]["load_bus"])
                                        input_ac_check["load"]["$(orig_f)"]["load_bus"] = deepcopy(input_dict["switch"]["$(switch_f["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["load"]["$(orig_f)"]["load_bus"]],"\n")
                                    elseif aux_f == "convdc"
                                        delete!(input_ac_check["bus"],input_ac_check["convdc"]["$(orig_f)"]["busac_i"])
                                        input_ac_check["convdc"]["$(orig_f)"]["busac_i"] = deepcopy(input_dict["switch"]["$(switch_f["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["convdc"]["$(orig_f)"]["busac_i"]],"\n")
                                    elseif aux_f == "branch"
                                        if input_ac_check["branch"]["$(orig_f)"]["f_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["f_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
                                            delete!(input_ac_check["bus"],input_ac_check["branch"]["$(orig_f)"]["f_bus"])
                                            input_ac_check["branch"]["$(orig_f)"]["f_bus"] = deepcopy(input_ac_check["switch"]["$(switch_f["index"])"]["t_bus"])
                                            print([l,aux_t,orig_t,input_ac_check["branch"]["$(orig_t)"]["f_bus"]],"\n")
                                        elseif input_ac_check["branch"]["$(orig_f)"]["t_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["f_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
                                            delete!(input_ac_check["bus"],input_ac_check["branch"]["$(orig_f)"]["t_bus"])
                                            input_ac_check["branch"]["$(orig_f)"]["t_bus"] = deepcopy(input_ac_check["switch"]["$(switch_f["index"])"]["t_bus"])
                                            print([l,aux_t,orig_t,input_ac_check["branch"]["$(orig_f)"]["t_bus"]],"\n")
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        input_ac_check["switch"] = Dict{String,Any}()
        input_ac_check["switch_couples"] = Dict{String,Any}()
end

function prepare_AC_feasibility_check(result_dict, input_dict, input_ac_check, switch_couples, extremes_dict,input_base)
    orig_buses = length(input_base["bus"]) # original bus length
    for (sw_id,sw) in input_ac_check["switch"]
        if !haskey(sw,"auxiliary")
            println("SWITCH $sw_id, BUS $(sw["t_bus"])")
            if result_dict["solution"]["switch"][sw_id]["status"] >= 0.9
                println("Switch $sw_id is closed, Connecting everything back, no busbar splitting on bus $(sw["bus_split"])")
                #delete!(input_ac_check["bus"],"$(input_ac_check["switch"][sw_id]["t_bus"])")
                for l in keys(switch_couples)
                    if switch_couples[l]["bus_split"] == sw["bus_split"]
                        println("SWITCH COUPLE IS $l")
                        println("Starting from switch $(switch_couples[l]["f_sw"]), with t_bus $(input_ac_check["switch"]["$(switch_couples[l]["f_sw"])"]["t_bus"])")
                        println("Then switch $(switch_couples[l]["t_sw"]), with t_bus $(input_ac_check["switch"]["$(switch_couples[l]["t_sw"])"]["t_bus"])")
                        if input_dict["switch"]["$(switch_couples[l]["f_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                            println("WE GO WITH SWITCH $(switch_couples[l]["f_sw"])")
                            aux =  deepcopy(input_dict["switch"]["$(switch_couples[l]["f_sw"])"]["auxiliary"])
                            orig = deepcopy(input_dict["switch"]["$(switch_couples[l]["f_sw"])"]["original"])
                            if aux == "gen"
                                input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "load"
                                input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "convdc"
                                input_ac_check["convdc"]["$(orig)"]["busac_i"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "branch"                
                                if input_dict["branch"]["$(orig)"]["f_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["f_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                    input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif input_dict["branch"]["$(orig)"]["t_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["f_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                    input_ac_check["branch"]["$(orig)"]["t_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                end
                            end
                        elseif input_dict["switch"]["$(switch_couples[l]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                            println("WE GO WITH SWITCH $(switch_couples[l]["t_sw"])")
                            println("----------------------------")
                            aux =  deepcopy(input_dict["switch"]["$(switch_couples[l]["t_sw"])"]["auxiliary"])
                            orig = deepcopy(input_dict["switch"]["$(switch_couples[l]["t_sw"])"]["original"])
                            if aux == "gen"
                                input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "load"
                                input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "convdc"
                                input_ac_check["convdc"]["$(orig)"]["busac_i"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "branch"                
                                if input_dict["branch"]["$(orig)"]["f_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                    input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif input_dict["branch"]["$(orig)"]["t_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                    input_ac_check["branch"]["$(orig)"]["t_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                end
                            end
                        end
                    end
                end
            elseif result_dict["solution"]["switch"][sw_id]["status"] <= 0.1
                println("Switch $sw_id is open, busbar splitting on bus $(sw["bus_split"])")
                delete!(input_ac_check["switch"],sw_id)
                for l in keys(switch_couples)
                    if switch_couples[l]["bus_split"] == sw["bus_split"]
                        println("SWITCH COUPLE IS $l")
                        println("Starting from switch $(switch_couples[l]["f_sw"]), with t_bus $(input_ac_check["switch"]["$(switch_couples[l]["f_sw"])"]["t_bus"])")
                        println("Then switch $(switch_couples[l]["t_sw"]), with t_bus $(input_ac_check["switch"]["$(switch_couples[l]["t_sw"])"]["t_bus"])")
                        if result_dict["solution"]["switch"]["$(switch_couples["$l"]["switch_split"])"]["status"] >= 0.9
                            aux =  deepcopy(input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["auxiliary"])
                            orig = deepcopy(input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["original"])
                            if aux == "gen"
                                input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "load"
                                input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "convdc"
                                input_ac_check["convdc"]["$(orig)"]["busac_i"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "branch"                
                                if input_dict["branch"]["$(orig)"]["f_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                    input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif input_dict["branch"]["$(orig)"]["t_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                    input_ac_check["branch"]["$(orig)"]["t_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                end
                            end
                        elseif result_dict["solution"]["switch"]["$(switch_couples["$l"]["switch_split"])"]["status"] <= 0.1
                            switch_t = deepcopy(input_ac_check["switch"]["$(switch_couples["$l"]["t_sw"])"]) 
                            aux_t = switch_t["auxiliary"]
                            orig_t = switch_t["original"]
                            print([l,aux_t,orig_t],"\n")
                            if result_dict["solution"]["switch"]["$(switch_t["index"])"]["status"] <= 0.1
                                delete!(input_ac_check["switch"],"$(switch_t["index"])")
                            elseif result_dict["solution"]["switch"]["$(switch_t["index"])"]["status"] >= 0.9
                                if aux_t == "gen"
                                    input_ac_check["gen"]["$(orig_t)"]["gen_bus"] = deepcopy(input_ac_check["switch"]["$(switch_t["index"])"]["t_bus"]) # here it needs to be the bus of the switch
                                    print([l,aux_t,orig_t,input_ac_check["gen"]["$(orig_t)"]["gen_bus"]],"\n")
                                elseif aux_t == "load"
                                    input_ac_check["load"]["$(orig_t)"]["load_bus"] = deepcopy(input_ac_check["switch"]["$(switch_t["index"])"]["t_bus"])
                                    print([l,aux_t,orig_t,input_ac_check["load"]["$(orig_t)"]["load_bus"]],"\n")
                                elseif aux_t == "convdc"
                                    input_ac_check["convdc"]["$(orig_t)"]["busac_i"] = deepcopy(input_ac_check["switch"]["$(switch_t["index"])"]["t_bus"])
                                    print([l,aux_t,orig_t,input_ac_check["convdc"]["$(orig_t)"]["busac_i"]],"\n")
                                elseif aux_t == "branch" 
                                    if input_dict["branch"]["$(orig_t)"]["f_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
                                        input_ac_check["branch"]["$(orig_t)"]["f_bus"] = deepcopy(input_ac_check["switch"]["$(switch_t["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["branch"]["$(orig_t)"]["f_bus"]],"\n")
                                    elseif input_dict["branch"]["$(orig_t)"]["t_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["t_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
                                        input_ac_check["branch"]["$(orig_t)"]["t_bus"] = deepcopy(input_ac_check["switch"]["$(switch_t["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["branch"]["$(orig_t)"]["t_bus"]],"\n")
                                    end
                                end
                            end
                        
                            switch_f = deepcopy(input_ac_check["switch"]["$(switch_couples["$l"]["f_sw"])"]) 
                            aux_f = switch_f["auxiliary"]
                            orig_f = switch_f["original"]
                            print([l,aux_f,orig_f],"\n")
                            if result_dict["solution"]["switch"]["$(switch_f["index"])"]["status"] <= 0.1
                                delete!(input_ac_check["switch"],"$(switch_t["index"])")
                            else
                                if aux_f == "gen"
                                    input_ac_check["gen"]["$(orig_f)"]["gen_bus"] = deepcopy(input_ac_check["switch"]["$(switch_f["index"])"]["t_bus"])
                                    print([l,aux_t,orig_t,input_ac_check["gen"]["$(orig_t)"]["gen_bus"]],"\n")
                                elseif aux_f == "load"
                                    input_ac_check["load"]["$(orig_f)"]["load_bus"] = deepcopy(input_ac_check["switch"]["$(switch_f["index"])"]["t_bus"])
                                    print([l,aux_t,orig_t,input_ac_check["load"]["$(orig_t)"]["load_bus"]],"\n")
                                elseif aux_f == "convdc"
                                    input_ac_check["convdc"]["$(orig_f)"]["busac_i"] = deepcopy(input_ac_check["switch"]["$(switch_f["index"])"]["t_bus"])
                                    print([l,aux_t,orig_t,input_ac_check["convdc"]["$(orig_t)"]["busac_i"]],"\n")
                                elseif aux_f == "branch"
                                    if input_dict["branch"]["$(orig_f)"]["f_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["f_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
                                        input_ac_check["branch"]["$(orig_f)"]["f_bus"] = deepcopy(input_ac_check["switch"]["$(switch_f["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["branch"]["$(orig_t)"]["f_bus"]],"\n")
                                    elseif input_dict["branch"]["$(orig_f)"]["t_bus"] > orig_buses && input_dict["switch"]["$(switch_couples["$l"]["f_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
                                        input_ac_check["branch"]["$(orig_f)"]["t_bus"] = deepcopy(input_ac_check["switch"]["$(switch_f["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["branch"]["$(orig_t)"]["t_bus"]],"\n")
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    delete!(input_ac_check,"switch")
    input_ac_check["switch"] = Dict()
end

###############################################################################
result_bs = Dict{String,Any}()
results_ac_check = Dict{String,Any}()
results_lpac_check = Dict{String,Any}()
split_one_bus_per_time(test_case,result_bs,results_ac_check,results_lpac_check)

objectives_bs = [result_bs["$b_id"]["objective"] for (b_id,b) in test_case["bus"] if result_bs["$b_id"]["termination_status"] == JuMP.OPTIMAL]
objectives_ac = [results_ac_check["$b_id"]["objective"] for (b_id,b) in test_case["bus"] if result_bs["$b_id"]["termination_status"] == JuMP.OPTIMAL]
objectives_lpac = [results_lpac_check["$b_id"]["objective"] for (b_id,b) in test_case["bus"] if result_bs["$b_id"]["termination_status"] == JuMP.OPTIMAL]

plot(objectives_bs)
plot!(objectives_ac)


result_bs_N_minus_1 = Dict{String,Any}()
results_ac_check_N_minus_1 = Dict{String,Any}()
results_lpac_check_N_minus_1 = Dict{String,Any}()

@time split_one_bus_per_time_with_N_1_check(test_case,result_bs_N_minus_1,results_ac_check_N_minus_1,results_lpac_check_N_minus_1)










result_bs_tree = Dict{String,Any}()
results_ac_check_tree = Dict{String,Any}()
results_lpac_check_tree = Dict{String,Any}()
busbar_splitting_tree(test_case,result_bs_tree,results_ac_check_tree,result_opf_lpac["objective"])


duals_plot = []
obj_plot = []
obj_plot_percentage = []
obj_plot_b_id = []

for (b_id,b) in test_case["bus"]
    if result_bs["$b_id"]["termination_status"] == JuMP.OPTIMAL
        push!(duals_plot,result_opf["solution"]["bus"][b_id]["lam_kcl_r"])
        push!(obj_plot,result_opf["objective"]-results_ac_check["$b_id"]["objective"])
        push!(obj_plot_percentage,((result_opf["objective"]-results_ac_check["$b_id"]["objective"])/result_opf["objective"])*100)
        push!(obj_plot_b_id,[b_id,result_opf["objective"]-results_ac_check["$b_id"]["objective"]])
    end
end

scatter(obj_plot_percentage,duals_plot,xlabel = "Benefit of busbar splitting compared to OPF [\$]",ylabel = "Duals [\$/MWh]",legend = false,grid = :none,ylims = (-5700,-2000),title = "Splitting one bus per time, case 30")
x_vertical_line = 0.0*ones(length(duals_plot))
vertical_line = -200*collect(1:length(duals_plot))
plot!(x_vertical_line,vertical_line,linecolor = :red)
savefig(joinpath(results_folder,"case_30/split_one_bus_per_time.pdf"))
savefig(joinpath(results_folder,"case_30/split_one_bus_per_time.svg"))

findmax(obj_plot)
sorted_bs = sort(obj_plot_b_id, by = x -> x[2], rev = true)

scatter(obj_plot_percentage,duals_plot,xlabel = "Benefit of busbar splitting compared to OPF [%]",ylabel = "Duals [\$/MWh]",legend = false,grid = :none,ylims = (-3700,-2000),title = "Splitting one bus per time")
x_vertical_line = 0.0*ones(length(duals_plot))
vertical_line = -200*collect(1:length(duals_plot))
plot!(x_vertical_line,vertical_line,linecolor = :red)



###########
results_folder = "/Users/giacomobastianel/Library/CloudStorage/OneDrive-KULeuven/Busbar_splitting_metrics_results"

bs_json = JSON.json(result_bs)
open(joinpath(results_folder,"case_30/split_one_bus_per_time.json"),"w") do f 
    write(f, bs_json) 
end
ac_check_json = JSON.json(results_ac_check)
open(joinpath(results_folder,"case_30/ac_check_split_one_bus_per_time.json"),"w") do f 
    write(f, ac_check_json) 
end
lpac_check_json = JSON.json(results_lpac_check)
open(joinpath(results_folder,"case_30/lpac_check_split_one_bus_per_time.json"),"w") do f 
    write(f, lpac_check_json) 
end

bs_json = JSON.json(result_bs_tree)
open(joinpath(results_folder,"case_30/split_tree.json"),"w") do f 
    write(f, bs_json) 
end
ac_check_json = JSON.json(results_ac_check_tree)
open(joinpath(results_folder,"case_30/ac_check_split_tree.json"),"w") do f 
    write(f, ac_check_json) 
end
lpac_check_json = JSON.json(results_lpac_check_tree)
open(joinpath(results_folder,"case_30/lpac_check_split_tree.json"),"w") do f 
    write(f, lpac_check_json) 
end

optimal_results = Dict{String,Any}()
optimal_results_ac_check = Dict{String,Any}()
optimal_results_lpac_check = Dict{String,Any}()
buses = []
for (b_id, result) in result_bs
    if result["termination_status"] == JuMP.OPTIMAL
        optimal_results[b_id] = result
        optimal_results_ac_check[b_id] = results_ac_check[b_id]
        optimal_results_lpac_check[b_id] = results_lpac_check[b_id]
        push!(buses,parse(Int64,b_id))
    end
end

obj_bs = [optimal_results["$b_id"]["objective"] for (b_id,b) in test_case["bus"] if result_bs["$b_id"]["termination_status"] == JuMP.OPTIMAL]
obj_bs_lpac_check = [optimal_results["$b_id"]["objective"] for (b_id,b) in test_case["bus"] if result_bs["$b_id"]["termination_status"] == JuMP.OPTIMAL]
obj = result_opf_lpac["objective"].*ones(length(obj_bs))

obj_bs_ac_check = [optimal_results_ac_check["$b_id"]["objective"] for (b_id,b) in test_case["bus"] if result_bs["$b_id"]["termination_status"] == JuMP.OPTIMAL]
obj_ac = result_opf["objective"].*ones(length(obj_bs))
