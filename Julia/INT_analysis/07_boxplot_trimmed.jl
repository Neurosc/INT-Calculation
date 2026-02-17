using Plots, StatsPlots, HypothesisTests, Statistics

# REMOVE OUTLIERS (tau > 30s)
# Some colors and fonts to use later on 
pltcolors= cgrad(:matter, 5, categorical = true)



# Function to remove outliers from subject means
function get_subject_mean_tau_trimmed(filepath, max_tau=30.0)
    data = load(filepath)
    scans = data["scans"]
    
    if isempty(scans)
        return NaN
    end
    
    all_tau_values = []
    for scan in scans
        tau_values = scan.acw_results[1]
        # Remove outliers before averaging
        tau_trimmed = tau_values[tau_values .<= max_tau]
        push!(all_tau_values, tau_trimmed)
    end
    
    if isempty(all_tau_values) || all(isempty.(all_tau_values))
        return NaN
    end
    
    # Flatten and average
    all_tau_flat = vcat(all_tau_values...)
    return mean(all_tau_flat)
end

# Recalculate for all subjects with trimming
control_self_trimmed = Float64[get_subject_mean_tau_trimmed(f, 30.0) for f in control_self_files]
control_nonself_trimmed = Float64[get_subject_mean_tau_trimmed(f, 30.0) for f in control_nonself_files]
med_self_trimmed = Float64[get_subject_mean_tau_trimmed(f, 30.0) for f in med_self_files]
med_nonself_trimmed = Float64[get_subject_mean_tau_trimmed(f, 30.0) for f in med_nonself_files]

# Clean NaNs
control_self_trimmed = filter(!isnan, control_self_trimmed)
control_nonself_trimmed = filter(!isnan, control_nonself_trimmed)
med_self_trimmed = filter(!isnan, med_self_trimmed)
med_nonself_trimmed = filter(!isnan, med_nonself_trimmed)

# Print summary
println("=== TRIMMED DATA (tau ≤ 30s) ===\n")
println("Controls Self: $(round(mean(control_self_trimmed), digits=2)) ± $(round(std(control_self_trimmed), digits=2))")
println("Controls Nonself: $(round(mean(control_nonself_trimmed), digits=2)) ± $(round(std(control_nonself_trimmed), digits=2))")
println("Meditators Self: $(round(mean(med_self_trimmed), digits=2)) ± $(round(std(med_self_trimmed), digits=2))")
println("Meditators Nonself: $(round(mean(med_nonself_trimmed), digits=2)) ± $(round(std(med_nonself_trimmed), digits=2))")

# STATISTICAL TESTS ON TRIMMED DATA

# Within-group
control_diff_trimmed = control_self_trimmed .- control_nonself_trimmed
med_diff_trimmed = med_self_trimmed .- med_nonself_trimmed

control_paired_trimmed = SignedRankTest(control_self_trimmed, control_nonself_trimmed)
med_paired_trimmed = SignedRankTest(med_self_trimmed, med_nonself_trimmed)

# Between-group
self_between_trimmed = MannWhitneyUTest(control_self_trimmed, med_self_trimmed)
nonself_between_trimmed = MannWhitneyUTest(control_nonself_trimmed, med_nonself_trimmed)

# Get p-values
p_control_paired_trimmed = pvalue(control_paired_trimmed)
p_med_paired_trimmed = pvalue(med_paired_trimmed)
p_self_between_trimmed = pvalue(self_between_trimmed)
p_nonself_between_trimmed = pvalue(nonself_between_trimmed)

# BOXPLOT WITH TRIMMED DATA

# Prepare data
control_self_labels_trimmed = fill("Control\nSelf", length(control_self_trimmed))
control_nonself_labels_trimmed = fill("Control\nNonself", length(control_nonself_trimmed))
med_self_labels_trimmed = fill("Meditator\nSelf", length(med_self_trimmed))
med_nonself_labels_trimmed = fill("Meditator\nNonself", length(med_nonself_trimmed))

all_data_trimmed = vcat(control_self_trimmed, control_nonself_trimmed,
                        med_self_trimmed, med_nonself_trimmed)
all_labels_trimmed = vcat(control_self_labels_trimmed, control_nonself_labels_trimmed,
                          med_self_labels_trimmed, med_nonself_labels_trimmed)

# Create boxplot
p = boxplot(all_labels_trimmed, all_data_trimmed,
        ylabel="Intrinsic Timescale (seconds)",
        title="Intrinsic Neural Timescales (Outliers Removed: τ > 30s)",
        legend=false,
        fillalpha=0.7,
        color= pltcolors[3],
        linewidth=2,
        size=(900, 700),
        xrotation=0,
        ylims=(0, maximum(all_data_trimmed) * 1.15))

# Add individual points
dotplot!(all_labels_trimmed, all_data_trimmed,
         marker=(:circle, 4, 0.4),
         color=:black)

# Add significance annotations
y_max = maximum(all_data_trimmed)

# Control: Self vs Nonself
annotate!(1.5, y_max * 1.05,
         text(p_control_paired_trimmed < 0.001 ? "***" : p_control_paired_trimmed < 0.01 ? "**" : p_control_paired_trimmed < 0.05 ? "*" : "ns",
              10, :center))
plot!([1, 2], [y_max * 1.02, y_max * 1.02], color=:black, linewidth=1.5)

# Meditator: Self vs Nonself
annotate!(3.5, y_max * 1.05,
         text(p_med_paired_trimmed < 0.001 ? "***" : p_med_paired_trimmed < 0.01 ? "**" : p_med_paired_trimmed < 0.05 ? "*" : "ns",
              10, :center))
plot!([3, 4], [y_max * 1.02, y_max * 1.02], color=:black, linewidth=1.5)

# Self: Control vs Meditator
annotate!(2, y_max * 1.12,
         text(p_self_between_trimmed < 0.001 ? "***" : p_self_between_trimmed < 0.01 ? "**" : p_self_between_trimmed < 0.05 ? "*" : "ns",
              10, :center))
plot!([1, 3], [y_max * 1.09, y_max * 1.09], color=:gray, linewidth=1.5)

# Nonself: Control vs Meditator
annotate!(3, y_max * 1.12,
         text(p_nonself_between_trimmed < 0.001 ? "***" : p_nonself_between_trimmed < 0.01 ? "**" : p_nonself_between_trimmed < 0.05 ? "*" : "ns",
              10, :center))
plot!([2, 4], [y_max * 1.09, y_max * 1.09], color=:gray, linewidth=1.5)

display(p)
savefig("timescales_boxplot_trimmed.png")

# PRINT STATISTICS

println("\n=== STATISTICAL TESTS (TRIMMED) ===")
println("\nWithin-Group:")
println("  Controls (Self vs Nonself): p = ", round(p_control_paired_trimmed, digits=4))
println("  Meditators (Self vs Nonself): p = ", round(p_med_paired_trimmed, digits=4))

println("\nBetween-Group:")
println("  Self ROIs: p = ", round(p_self_between_trimmed, digits=4))
println("  Nonself ROIs: p = ", round(p_nonself_between_trimmed, digits=4))