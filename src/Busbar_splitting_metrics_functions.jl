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
            results_dict_ac_check["$b_id"] = _PMACDC.run_acdcopf(test_case_bs_check,ACPPowerModel,ipopt; setting = s)
            results_dict_lpac_check["$b_id"] = _PMACDC.run_acdcopf(test_case_bs_check,LPACCPowerModel,gurobi; setting = s)
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