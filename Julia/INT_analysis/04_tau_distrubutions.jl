# This file examines the distribution of tau values pooled across all subjects and scans. 
# Results are visualized as percentage histograms to compare self and nonself ROI distributions.


using Plots, Statistics

# COLLECT ALL TAU VALUES ACROSS ALL SUBJECTS

# Function to get all tau values from all scans for a subject
function get_all_tau_values(filepath)
    data = load(filepath)
    scans = data["scans"]
    
    all_tau = Float64[]
    for scan in scans
        if !isempty(scan.acw_results)
            append!(all_tau, scan.acw_results[1])
        end
    end
    return all_tau
end

# Collect all self tau values (across all controls and meditators)
all_self_tau = Float64[]
for file in vcat(control_self_files, med_self_files)
    append!(all_self_tau, get_all_tau_values(file))
end

# Collect all nonself tau values
all_nonself_tau = Float64[]
for file in vcat(control_nonself_files, med_nonself_files)
    append!(all_nonself_tau, get_all_tau_values(file))
end

# SUMMARY STATISTICS

println("=== DISTRIBUTION COMPARISON ===\n")

println("Self ROIs (pooled across all subjects):")
println("  Total measurements: ", length(all_self_tau))
println("  Mean: ", round(mean(all_self_tau), digits=2))
println("  Median: ", round(median(all_self_tau), digits=2))
println("  Max: ", round(maximum(all_self_tau), digits=2))
println("  ROIs > 30s: ", sum(all_self_tau .> 30), " (", round(sum(all_self_tau .> 30)/length(all_self_tau)*100, digits=2), "%)")
println("  ROIs > 40s: ", sum(all_self_tau .> 40), " (", round(sum(all_self_tau .> 40)/length(all_self_tau)*100, digits=2), "%)")

println("\nNonself ROIs (pooled across all subjects):")
println("  Total measurements: ", length(all_nonself_tau))
println("  Mean: ", round(mean(all_nonself_tau), digits=2))
println("  Median: ", round(median(all_nonself_tau), digits=2))
println("  Max: ", round(maximum(all_nonself_tau), digits=2))
println("  ROIs > 30s: ", sum(all_nonself_tau .> 30), " (", round(sum(all_nonself_tau .> 30)/length(all_nonself_tau)*100, digits=2), "%)")
println("  ROIs > 40s: ", sum(all_nonself_tau .> 40), " (", round(sum(all_nonself_tau .> 40)/length(all_nonself_tau)*100, digits=2), "%)")
