# This file contains all reusable functions for loading and processing ACW results.

# Functions include:
#   1. get_subject_tau(filepath; average_rois=true, use_median=false)
#      - Returns tau values for a subject
#      - average_rois=true :(mean/median across all ROIs)
#      - average_rois=false :returns vector of tau per ROI (e.g., 37 or 327 values)
#      - use_median=false   :uses mean
#      - use_median=true    :uses median

#   2. get_subject_tau_bootstrap(filepath, n_bootstrap=1000; use_median=false)
#      - For NONSELF ROIs only (327 ROIs)
#      - Randomly samples 37 ROIs, repeats n_bootstrap times
#      - Returns average of bootstrap samples
#      - use_median=false :uses mean
#      - use_median=true  :uses median

using JLD2, Statistics, StatsBase, Random, StatsPlots, Plots, ColorSchemes

function get_subject_tau(filepath; average_rois=true, use_median=false)
    data = load(filepath)
    scans = data["scans"]

    if isempty(scans)
        return average_rois ? NaN : nothing
    end

    all_tau_values = Vector{Vector{Float64}}()
    for scan in scans
        push!(all_tau_values, scan.acw_results[1])
    end

    if isempty(all_tau_values)
        return average_rois ? NaN : nothing
    end

    # Stack: scans × ROIs
    tau_matrix = reduce(hcat, all_tau_values)'

    # Average or median across scans
    summary_per_roi = use_median ? vec(median(tau_matrix, dims=1)) : vec(mean(tau_matrix, dims=1))

    # Return either one number or all ROIs
    if average_rois
        return use_median ? median(summary_per_roi) : mean(summary_per_roi)
    else
        return summary_per_roi
    end
end


function get_subject_tau_bootstrap(filepath, n_bootstrap=1000; use_median=false)
    data = load(filepath)
    scans = data["scans"]

    if isempty(scans)
        return NaN
    end

    all_tau_values = Vector{Vector{Float64}}()
    for scan in scans
        push!(all_tau_values, scan.acw_results[1])
    end

    if isempty(all_tau_values)
        return NaN
    end

    # Stack: scans × ROIs
    tau_matrix = reduce(hcat, all_tau_values)'

    # Average or median across scans
    summary_per_roi = use_median ? vec(median(tau_matrix, dims=1)) : vec(mean(tau_matrix, dims=1))

    n_rois = length(summary_per_roi)
    bootstrap_values = Float64[]

    for i in 1:n_bootstrap
        sampled_rois = sample(1:n_rois, 37, replace=false)
        val = use_median ? median(summary_per_roi[sampled_rois]) : mean(summary_per_roi[sampled_rois])
        push!(bootstrap_values, val)
    end

    return use_median ? median(bootstrap_values) : mean(bootstrap_values)
end