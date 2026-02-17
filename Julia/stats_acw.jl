using JLD2, Statistics, StatsBase, Random, StatsPlots, Plots, ColorSchemes


function get_subject_mean_tau(filepath)
    data = load(filepath)
    scans = data["scans"]

    if isempty(scans)
        return NaN
    end

    all_tau_values = Vector{Vector{Float64}}()

    for scan in scans   
        push!(all_tau_values, scan.acw_results[1]) # auc didnt work.
    end

    if isempty(all_tau_values)
        return NaN
    end

    # Stack: scans × ROIs
    tau_matrix = reduce(hcat, all_tau_values)' # hcat concatenates them, then we transpose so it s scansxROIs
   
    # Average across scans, then across ROIs
    mean_per_roi = vec(mean(tau_matrix, dims=1)) #vec makes it 37-element

    return mean(mean_per_roi)
end


function get_subject_mean_tau_bootstrap(filepath, n_bootstrap=1000)
        data = load(filepath)
    scans = data["scans"]

    if isempty(scans)
        return NaN
    end

    all_tau_values = Vector{Vector{Float64}}()

    for scan in scans   
        push!(all_tau_values, scan.acw_results[1]) # auc didnt work.
    end

    if isempty(all_tau_values)
        return NaN
    end

    # Stack: scans × ROIs
    tau_matrix = reduce(hcat, all_tau_values)' # hcat concatenates them, then we transpose so it s scansxROIs
   
    mean_per_roi = vec(mean(tau_matrix, dims=1))

    n_rois = length(mean_per_roi)
    bootstrap_means = Float64[]

    for i in 1:n_bootstrap
        sampled_rois= sample(1:n_rois, 37, replace=false) # pick any number in between 1:327, pick 37 of them, dont use the same number twice
        push!(bootstrap_means, mean(mean_per_roi[sampled_rois])) # then take those indice values from mean_per_roi, and put it into bootstrap_means
    end

    return mean(bootstrap_means)
end

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

# MEDITATORS
med_files = readdir(meditator_path, join=true)
med_self_files = filter(f -> contains(f, "extracted_self") && endswith(f, ".jld2"), med_files)
med_nonself_files = filter(f -> contains(f, "extracted_nonself") && endswith(f, ".jld2"), med_files)

med_self = Float64[get_subject_mean_tau(f) for f in med_self_files]
med_nonself = Float64[get_subject_mean_tau_bootstrap(f) for f in med_nonself_files]

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

# STEP 2: CREATE BOXPLOT WITH STATISTICS

# Some colors and fonts to use later on 
pltcolors= cgrad(:matter, 5, categorical = true)
default(fontfamily="Computer Modern")  # LaTeX-style font


# Prepare data
control_self_labels = fill("Control\nSelf", length(control_self_clean))
control_nonself_labels = fill("Control\nNonself", length(control_nonself_clean))
med_self_labels = fill("Meditator\nSelf", length(med_self_clean))
med_nonself_labels = fill("Meditator\nNonself", length(med_nonself_clean))

all_data = vcat(control_self_clean, control_nonself_clean, 
                med_self_clean, med_nonself_clean)
all_labels = vcat(control_self_labels, control_nonself_labels,
                  med_self_labels, med_nonself_labels)
# Create boxplot
p = boxplot(all_labels, all_data,
        ylabel="Intrinsic Timescale (seconds)",
        title="Intrinsic Neural Timescales: Self vs Nonself ROIs",
        legend=false,
        fillalpha=0.7,
        linewidth=2,
        size=(900, 700),
        xrotation=0,
        ylims=(0, maximum(all_data) * 1.15),
        color= pltcolors[3])  # Extra space for annotations

# Add individual points
dotplot!(all_labels, all_data, 
         marker=(:circle, 4, 0.4),
         color=pltcolors[1])



display(p)
savefig("timescales_boxplot_with_stats.png")


# STEP 3: PRINT STATISTICS SUMMARY

println("\n=== STATISTICAL SUMMARY ===")
println("\nWithin-Group (Wilcoxon signed-rank):")
println("  Controls (Self vs Nonself): p = ", round(p_control_paired, digits=4))
println("  Meditators (Self vs Nonself): p = ", round(p_med_paired, digits=4))

println("\nBetween-Group (Mann-Whitney U):")
println("  Self ROIs (C vs M): p = ", round(p_self_between, digits=4))
println("  Nonself ROIs (C vs M): p = ", round(p_nonself_between, digits=4))

println("\n* p < 0.05, ** p < 0.01, *** p < 0.001")