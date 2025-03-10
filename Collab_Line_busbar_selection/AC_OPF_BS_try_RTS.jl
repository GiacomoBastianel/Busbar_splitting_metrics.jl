using PowerModels; const _PM = PowerModels
using PowerModelsACDC; const _PMACDC = PowerModelsACDC
using PowerModelsTopologicalActionsII; const _PMTP = PowerModelsTopologicalActionsII
using Gurobi, JuMP, HiGHS, DataFrames, CSV, Feather, JSON, Ipopt, Plots, StatsBase, Juniper, Statistics, PowerPlots
#import HSL_jll


gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"QCPDual" => 1, "time_limit" => 600,"MIPgap" =>1e-3, "BarHomogeneous" => 1,"BarQCPConvTol"=>1e-4)#r, "ScaleFlag"=>2, "NumericFocus"=>2,,"BarQCPConvTol"=>1e-3,) 
gurobi_lpac = JuMP.optimizer_with_attributes(Gurobi.Optimizer)

ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 2)#,"linear_solver" => "ma97")
juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt, "mip_solver" => gurobi, "time_limit" => 1200, "allow_almost_solved_integral" => false)
s = Dict("output" => Dict("duals" => false))
s_dual = Dict("output" => Dict("duals" => true))

##################################################################
## Processing input data (RTS-GMLC)
input_folder = dirname(@__DIR__)
test_case_file = joinpath(input_folder,"test_cases/RTS_GMLC.m")
test_case = _PM.parse_file(test_case_file)
test_case_lpac = _PM.parse_file(test_case_file)

result_opf = _PM.solve_opf(test_case,ACPPowerModel,ipopt,setting=s_dual)

duals = []
for (b_id,b) in test_case["bus"]
    push!(duals,[b_id,result_opf["solution"]["bus"][b_id]["lam_kcl_r"]])
end
sort!(duals, by = x -> x[2], rev = false)

result_opf_lpac = _PM.solve_opf(test_case_lpac,LPACCPowerModel,gurobi_lpac,setting=s)
result_opf_soc = _PM.solve_opf(test_case_lpac,SOCWRPowerModel,gurobi_lpac,setting=s)
result_opf_dc = _PM.solve_opf(test_case,DCPPowerModel,gurobi,setting=s)
result_opf_ac = _PM.solve_opf(test_case,ACPPowerModel,ipopt,setting=s)


#######################################################
result_bs = Dict{String,Any}()
results_ac_check = Dict{String,Any}()
results_lpac_check = Dict{String,Any}()
split_one_bus_per_time(test_case,result_bs,results_ac_check,results_lpac_check)
objectives_bs = [result_bs["$b_id"]["objective"] for (b_id,b) in test_case["bus"] if result_bs["$b_id"]["termination_status"] == JuMP.OPTIMAL]


duals_plot = []
obj_plot = []
obj_plot_percentage = []
obj_plot_b_id = []

for (b_id,b) in test_case["bus"]
    if result_bs["$b_id"]["termination_status"] == JuMP.OPTIMAL
        push!(duals_plot,result_opf["solution"]["bus"][b_id]["lam_kcl_r"])
        push!(obj_plot,result_opf["objective"]-results_ac_check["$b_id"]["objective"])
        push!(obj_plot_percentage,((result_opf["objective"]-results_ac_check["$b_id"]["objective"])/result_opf["objective"])*100)
        push!(obj_plot_b_id,[b_id,result_opf_lpac["objective"]-result_bs["$b_id"]["objective"]])
    end
end


scatter(obj_plot,duals_plot,xlabel = "Benefit of busbar splitting compared to OPF [\$]",ylabel = "Duals [\$/MWh]",legend = false,grid = :none,ylims = (-5700,-2000),title = "Splitting one bus per time, case RTS GMLC")
x_vertical_line = 0.0*ones(length(duals_plot))
vertical_line = -200*collect(1:length(duals_plot))
plot!(x_vertical_line,vertical_line,linecolor = :red)
savefig(joinpath(results_folder,"case_118/split_one_bus_per_time.pdf"))
savefig(joinpath(results_folder,"case_118/split_one_bus_per_time.svg"))
















buses = []
bs_opf_diff = []
dual_buses = []
for (b_id,b) in test_case_bs_try["bus"]
    push!(buses,parse(Int64,b_id))
    push!(dual_buses,result_opf["solution"]["bus"][b_id]["lam_kcl_r"])
    push!(bs_opf_diff,result_opf["objective"] - result[b_id]["objective"])
end
zero = 0.0*collect(1:length(bs_opf_diff))
vertical_line = -100*collect(1:length(bs_opf_diff))
scatter(bs_opf_diff/10^3,dual_buses,legend= :none,grid = :none,ylabel = "Bus' dual (lam_kcl_r) in OPF",xlabel = "Diff AC-OPF - AC-feasibility check [k€], > 0 = reduction in costs",ylims = (-5000,-3000),xlims = (-30,15),title = "RTS GMC test case, busbar splitting, no OTS")
annotate!([(bs_opf_diff[i]/10^3, dual_buses[i], text(buses[i], 9, :red, :right)) for i in 1:length(buses) if bs_opf_diff[i] > 0])
plot!(zero,vertical_line,legend= :none,color = :red)





