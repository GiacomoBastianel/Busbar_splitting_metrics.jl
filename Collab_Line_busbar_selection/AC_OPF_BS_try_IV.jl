using PowerModels; const _PM = PowerModels
using PowerModelsACDC; const _PMACDC = PowerModelsACDC
using PowerModelsTopologicalActionsII; const _PMTP = PowerModelsTopologicalActionsII
using Gurobi, JuMP, HiGHS, DataFrames, CSV, Feather, JSON, Ipopt, Plots, StatsBase, Juniper, Statistics, PowerPlots

mip_gap = 5e-4
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"QCPDual" => 1, "time_limit" => 600,"MIPgap" => mip_gap, "BarHomogeneous" => 1,"BarQCPConvTol"=>1e-4)#r, "ScaleFlag"=>2, "NumericFocus"=>2,,"BarQCPConvTol"=>1e-3,) 
#ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 2,"linear_solver" => "ma97")
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 2)
juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt, "mip_solver" => gurobi, "time_limit" => 600)
highs = JuMP.optimizer_with_attributes(HiGHS.Optimizer)#, "presolve" => 0, "presolve_dual" => 0, "simplex_dual_edge_weight_strategy" => 0, "simplex_primal_edge_weight_strategy" => 0, "simplex_price_strategy" => 0, "simplex_iteration_limit" => 100000, "simplex_perturbation" => 50.0, "simplex_dual_feasibility_tolerance" => 1.0e-9, "simplex_primal_feasibility_tolerance" => 1.0e-9, "simplex_objective_tolerance" => 1.0e-9, "simplex_infinite_bound" => 1.0e20, "simplex_infinite_cost" => 1.0e20)

##################################################################
## Processing input data
input_folder = @__DIR__

# Belgium grid without energy island
test_case_file = joinpath(dirname(input_folder),"test_cases/pglib_opf_case118_ieee.m")
test_case = _PM.parse_file(test_case_file)

s_dual = Dict("output" => Dict("duals" => true))
s = Dict("output" => Dict("duals" => false))

result_opf = _PM.solve_opf(test_case,ACPPowerModel,ipopt,setting=s)
result_opf_lpac = _PM.solve_opf(test_case,LPACCPowerModel,gurobi,setting=s)
#result_opf_soc = _PM.solve_opf(test_case,SOCWRPowerModel,gurobi,setting=s)


powerplot(test_case)

for (b_id,b) in test_case["bus"]
    println(b_id," ",result_opf["solution"]["bus"][b_id]["lam_kcl_r"])
end
duals = []
for (b_id,b) in test_case["bus"]
    push!(duals,[parse(Int64,b_id),abs(result_opf["solution"]["bus"][b_id]["lam_kcl_r"])])
end

sorted_duals = sort(duals, by = x -> x[2], rev = true)

#=
for (b_id,b) in test_case["branch"]
    println(b_id," ",result_opf["solution"]["branch"][b_id]["pf"]/b["rate_a"])
end

for (b_id,b) in test_case["branch"]
    println(b_id," f_bus ",b["f_bus"]," t_bus ",b["t_bus"]," ",abs(abs(result_opf["solution"]["bus"]["$(b["f_bus"])"]["lam_kcl_r"]) - abs(result_opf["solution"]["bus"]["$(b["t_bus"])"]["lam_kcl_r"])))
end
duals = []
for (b_id,b) in test_case["branch"]
    push!(duals,[parse(Int64,b_id),abs(abs(result_opf["solution"]["bus"]["$(b["f_bus"])"]["lam_kcl_r"]) - abs(result_opf["solution"]["bus"]["$(b["t_bus"])"]["lam_kcl_r"]))])
end
sorted_duals = sort(duals, by = x -> x[2], rev = true)
=#
test_case_bs = deepcopy(test_case)

# Selecting which busbars are split
splitted_bus_ac = 49
#splitted_bus_ac = [69,103,32,77,59] #-> simulation: 95s, N-1 check 22.28s (all the lines! + everything saved in Dict) -> < 2 mins vs 4+ mins comparable to Heidarifar et al. 2021

optimizer = gurobi
formulation = LPACCPowerModel

test_case_bs,  switches_couples_ac,  extremes_ZILs_ac  = _PMTP.AC_busbar_split_AC_grid(test_case_bs,splitted_bus_ac)
test_case_bs_check = deepcopy(test_case_bs)
test_case_bs_check_auxiliary = deepcopy(test_case_bs)
results_bs = _PMTP.run_acdcsw_AC_grid_big_M(test_case_bs,formulation,optimizer)

test_case_bs_sp = deepcopy(test_case_bs)
function prepare_starting_value_dict_lpac(result,grid)
    for (b_id,b) in grid["bus"]
        if haskey(result["solution"]["bus"],b_id)
            if abs(result["solution"]["bus"]["$b_id"]["va"]) < 10^(-4)
                b["va_starting_value"] = 0.0
            else
                b["va_starting_value"] = result["solution"]["bus"]["$b_id"]["va"]
            end
            if abs(result["solution"]["bus"]["$b_id"]["phi"]) < 10^(-4)
                b["phi_starting_value"] = 0.0
            else
                b["phi_starting_value"] = result["solution"]["bus"]["$b_id"]["phi"]
            end
        else
            b["va_starting_value"] = 0.0
            b["phi_starting_value"] = 0.1
        end
    end
    for (b_id,b) in grid["gen"]
        if abs(result["solution"]["gen"]["$b_id"]["pg"]) < 10^(-5)
            b["pg_starting_value"] = 0.0
        else
            b["pg_starting_value"] = result["solution"]["gen"]["$b_id"]["pg"]
        end
        if abs(result["solution"]["gen"]["$b_id"]["qg"]) < 10^(-5)
            b["qg_starting_value"] = 0.0
        else
            b["qg_starting_value"] = result["solution"]["gen"]["$b_id"]["qg"]
        end
    end
    for (sw_id,sw) in grid["switch"]
        if !haskey(sw,"auxiliary") # calling ZILs
            sw["starting_value"] = 1.0
        else
            if haskey(grid["switch_couples"],sw_id)
                grid["switch"]["$(grid["switch_couples"][sw_id]["f_sw"])"]["starting_value"] = 0.0
                grid["switch"]["$(grid["switch_couples"][sw_id]["t_sw"])"]["starting_value"] = 1.0
            end
        end
    end
end
prepare_starting_value_dict_lpac(result_opf_lpac,test_case_bs_sp)

results_bs_sp = _PMTP.run_acdcsw_AC_grid_big_M_sp(test_case_bs,formulation,optimizer)


function prepare_AC_grid_feasibility_check(result_dict, input_dict, input_ac_check, switch_couples, extremes_dict,input_base)
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
                        if input_ac_check["switch"]["$(switch_couples[l]["f_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                            println("WE GO WITH SWITCH $(switch_couples[l]["f_sw"])")
                            aux =  deepcopy(input_ac_check["switch"]["$(switch_couples[l]["f_sw"])"]["auxiliary"])
                            orig = deepcopy(input_ac_check["switch"]["$(switch_couples[l]["f_sw"])"]["original"])
                            if aux == "gen"
                                input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "load"
                                input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "convdc"
                                input_ac_check["convdc"]["$(orig)"]["busac_i"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "branch"                
                                if input_ac_check["branch"]["$(orig)"]["f_bus"] > orig_buses && input_ac_check["switch"]["$(switch_couples["$l"]["f_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                    input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif input_ac_check["branch"]["$(orig)"]["t_bus"] > orig_buses && input_ac_check["switch"]["$(switch_couples["$l"]["f_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                    input_ac_check["branch"]["$(orig)"]["t_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                end
                            end
                        elseif input_ac_check["switch"]["$(switch_couples[l]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                            println("WE GO WITH SWITCH $(switch_couples[l]["t_sw"])")
                            println("----------------------------")
                            aux =  deepcopy(input_ac_check["switch"]["$(switch_couples[l]["t_sw"])"]["auxiliary"])
                            orig = deepcopy(input_ac_check["switch"]["$(switch_couples[l]["t_sw"])"]["original"])
                            if aux == "gen"
                                input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "load"
                                input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "convdc"
                                input_ac_check["convdc"]["$(orig)"]["busac_i"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "branch"                
                                if input_ac_check["branch"]["$(orig)"]["f_bus"] > orig_buses && input_ac_check["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                    input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif input_ac_check["branch"]["$(orig)"]["t_bus"] > orig_buses && input_ac_check["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
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
                            aux =  deepcopy(input_ac_check["switch"]["$(switch_couples["$l"]["t_sw"])"]["auxiliary"])
                            orig = deepcopy(input_ac_check["switch"]["$(switch_couples["$l"]["t_sw"])"]["original"])
                            if aux == "gen"
                                input_ac_check["gen"]["$(orig)"]["gen_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "load"
                                input_ac_check["load"]["$(orig)"]["load_bus"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "convdc"
                                input_ac_check["convdc"]["$(orig)"]["busac_i"] = deepcopy(switch_couples[l]["bus_split"])
                            elseif aux == "branch"                
                                if input_ac_check["branch"]["$(orig)"]["f_bus"] > orig_buses && input_ac_check["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
                                    input_ac_check["branch"]["$(orig)"]["f_bus"] = deepcopy(switch_couples[l]["bus_split"])
                                elseif input_ac_check["branch"]["$(orig)"]["t_bus"] > orig_buses && input_ac_check["switch"]["$(switch_couples["$l"]["t_sw"])"]["t_bus"] == switch_couples[l]["bus_split"]
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
                                    if input_ac_check["branch"]["$(orig_t)"]["f_bus"] > orig_buses && input_ac_check["switch"]["$(switch_couples["$l"]["t_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
                                        input_ac_check["branch"]["$(orig_t)"]["f_bus"] = deepcopy(input_ac_check["switch"]["$(switch_t["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["branch"]["$(orig_t)"]["f_bus"]],"\n")
                                    elseif input_ac_check["branch"]["$(orig_t)"]["t_bus"] > orig_buses && input_ac_check["switch"]["$(switch_couples["$l"]["t_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
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
                                    if input_ac_check["branch"]["$(orig_f)"]["f_bus"] > orig_buses && input_ac_check["switch"]["$(switch_couples["$l"]["f_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
                                        input_ac_check["branch"]["$(orig_f)"]["f_bus"] = deepcopy(input_ac_check["switch"]["$(switch_f["index"])"]["t_bus"])
                                        print([l,aux_t,orig_t,input_ac_check["branch"]["$(orig_t)"]["f_bus"]],"\n")
                                    elseif input_ac_check["branch"]["$(orig_f)"]["t_bus"] > orig_buses && input_ac_check["switch"]["$(switch_couples["$l"]["f_sw"])"]["bus_split"] == switch_couples[l]["bus_split"]
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

test_case_bs_check = deepcopy(test_case_bs)
prepare_AC_feasibility_check(results_bs,test_case_bs_check_auxiliary,test_case_bs_check,switches_couples_ac,extremes_ZILs_ac,test_case)
result_opf_check = _PMACDC.run_acdcopf(test_case_bs_check,ACPPowerModel,ipopt; setting = s)

println("Benefits from splitting busbar 49 is $((result_opf["objective"] - result_opf_check["objective"])), $((result_opf["objective"]-result_opf_check["objective"])/(result_opf["objective"])*100) %, with MIP gap $(mip_gap*100) %")


function N_minus_1_checks(grid,results_dict)
    grid_base = deepcopy(grid)
    for (br_id,br) in grid["branch"]
        hourly_grid = deepcopy(grid_base)
        br["br_status"] = 0
        hourly_result = _PM.solve_opf(hourly_grid,ACPPowerModel,ipopt,setting=s)
        if length(hourly_result["solution"]) <= 5
            results_dict["$(br_id)"] = deepcopy(hourly_result)
        end
    end
    return results_dict
end
Results_N_minus_1 = Dict{String,Any}()
@time N_minus_1_checks(test_case_bs_check,Results_N_minus_1)


powerplot(test_case_bs_check)

##################################################################################
splitted_bus_ac = [69,103,32,77,59] #-> simulation: 95s, N-1 check 22.28s (all the lines! + everything saved in Dict) -> < 2 mins vs 4+ mins comparable to Heidarifar et al. 2021

optimizer = gurobi
formulation = LPACCPowerModel

test_case_bs = deepcopy(test_case)
test_case_bs,  switches_couples_ac,  extremes_ZILs_ac  = _PMTP.AC_busbar_split_AC_grid(test_case_bs,splitted_bus_ac)
test_case_bs_check = deepcopy(test_case_bs)
test_case_bs_check_auxiliary = deepcopy(test_case_bs)
results_bs = _PMTP.run_acdcsw_AC_grid_big_M(test_case_bs,formulation,optimizer)
results_bs_so = _PMTP.run_acdcsw_AC_grid_big_M_sp(test_case_bs,formulation,optimizer)


prepare_AC_grid_feasibility_check(results_bs,test_case_bs_check_auxiliary,test_case_bs_check,switches_couples_ac,extremes_ZILs_ac,test_case)
result_opf_check = _PMACDC.run_acdcopf(test_case_bs_check,ACPPowerModel,ipopt; setting = s)

println("Benefits from splitting busbar 69,103,32,77,59 is $((result_opf["objective"] - result_opf_check["objective"])), $((result_opf["objective"]-result_opf_check["objective"])/(result_opf["objective"])*100)%, with MIPgap $(mip_gap*100) %")

Results_N_minus_1 = Dict{String,Any}()
@time N_minus_1_checks(test_case_bs_check,Results_N_minus_1)


obj = [Results_N_minus_1["$br_id"]["objective"] for (br_id,br) in test_case_bs_check["branch"]]
termination_status = [Results_N_minus_1["$br_id"]["termination_status"] for (br_id,br) in test_case_bs_check["branch"]]

powerplot(test_case_bs_check)
