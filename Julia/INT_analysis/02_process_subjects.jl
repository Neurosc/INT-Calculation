# This file loads and processes ACW results for all control and meditator subjects. 

# Outputs are:
#### cleaned group-level mean tau values 
#### for self and nonself ROIs.


Random.seed!(123)

# Paths
control_path = "C:/Users/JLU-MBB/Desktop/Thesis/ACW_results/controls"
meditator_path = "C:/Users/JLU-MBB/Desktop/Thesis/ACW_results/meditators"

# CONTROLS
control_files = readdir(control_path, join=true)
control_self_files = filter(f -> contains(f, "extracted_self") && endswith(f, ".jld2"), control_files)
control_nonself_files = filter(f -> contains(f, "extracted_nonself") && endswith(f, ".jld2"), control_files)

control_self = Float64[get_subject_mean_tau(f) for f in control_self_files]
control_nonself = Float64[get_subject_mean_tau_bootstrap(f) for f in control_nonself_files]
#control_self = Float64[get_subject_median_tau(f) for f in control_self_files]
#control_nonself = Float64[get_subject_median_tau_bootstrap(f) for f in control_nonself_files]


# MEDITATORS
med_files = readdir(meditator_path, join=true)
med_self_files = filter(f -> contains(f, "extracted_self") && endswith(f, ".jld2"), med_files)
med_nonself_files = filter(f -> contains(f, "extracted_nonself") && endswith(f, ".jld2"), med_files)

med_self = Float64[get_subject_mean_tau(f) for f in med_self_files]
med_nonself = Float64[get_subject_mean_tau_bootstrap(f) for f in med_nonself_files]
#med_self = Float64[get_subject_median_tau(f) for f in med_self_files]
#med_nonself = Float64[get_subject_median_tau_bootstrap(f) for f in med_nonself_files]

# Clean NaNs
control_self_clean = filter(!isnan, control_self)
control_nonself_clean = filter(!isnan, control_nonself)
med_self_clean = filter(!isnan, med_self)
med_nonself_clean = filter(!isnan, med_nonself)

# Summary
println("=== SUMMARY ===")
println("Controls Self (N=$(length(control_self_clean))): $(round(mean(control_self_clean), digits=2)) ± $(round(std(control_self_clean), digits=2))")
println("Controls Nonself (N=$(length(control_nonself_clean))): $(round(mean(control_nonself_clean), digits=2)) ± $(round(std(control_nonself_clean), digits=2))")
println("Meditators Self (N=$(length(med_self_clean))): $(round(mean(med_self_clean), digits=2)) ± $(round(std(med_self_clean), digits=2))")
println("Meditators Nonself (N=$(length(med_nonself_clean))): $(round(mean(med_nonself_clean), digits=2)) ± $(round(std(med_nonself_clean), digits=2))")

using Plots, StatsPlots, HypothesisTests, Statistics

# Step 1: Statistics

#Within-group comparisons (Paired: Self vs Nonself)
control_diff = control_self_clean .- control_nonself_clean
med_diff = med_self_clean .- med_nonself_clean

control_paired = SignedRankTest(control_self_clean, control_nonself_clean)
med_paired = SignedRankTest(med_self_clean, med_nonself_clean)

# Between-group comparisons (Independent: Controls vs Meditators)
self_between = MannWhitneyUTest(control_self_clean, med_self_clean)
nonself_between = MannWhitneyUTest(control_nonself_clean, med_nonself_clean)

# Get p-values
p_control_paired = pvalue(control_paired)
p_med_paired = pvalue(med_paired)
p_self_between = pvalue(self_between)
p_nonself_between = pvalue(nonself_between)


# Combine self and nonself for each group
control_all = vcat(control_self_clean, control_nonself_clean)
meditator_all = vcat(med_self_clean, med_nonself_clean)

println("=== OVERALL MEDITATION EFFECT ===\n")
println("Controls (N=$(length(control_all))): $(round(mean(control_all), digits=2)) ± $(round(std(control_all), digits=2))")
println("Meditators (N=$(length(meditator_all))): $(round(mean(meditator_all), digits=2)) ± $(round(std(meditator_all), digits=2))")

# Mann-Whitney U test
overall_test = MannWhitneyUTest(control_all, meditator_all)
println("\nMann-Whitney U: p = ", round(pvalue(overall_test), digits=4))

# Effect size
pooled_std = sqrt((std(control_all)^2 + std(meditator_all)^2) / 2)
cohens_d = (mean(control_all) - mean(meditator_all)) / pooled_std
println("Cohen's d = ", round(cohens_d, digits=3))