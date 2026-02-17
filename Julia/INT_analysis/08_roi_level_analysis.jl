# PERMUTATION TEST: ALL 364 ROIs (MEDIAN)

n_permutations = 1000
Random.seed!(123)

perm_pvalues_med = Float64[]
perm_differences_med = Float64[]


for roi in 1:size(control_combined_matrix, 2)
    control_roi = control_combined_matrix[:, roi]
    med_roi = med_combined_matrix[:, roi]
    
    # Real difference (median)
    real_diff = median(control_roi) - median(med_roi)
    push!(perm_differences_med, real_diff)
    
    # Permutation test
    all_subjects = vcat(control_roi, med_roi)
    n_controls = length(control_roi)
    
    null_distribution = Float64[]
    for i in 1:n_permutations
        shuffled = shuffle(all_subjects)
        perm_diff = median(shuffled[1:n_controls]) - median(shuffled[n_controls+1:end])
        push!(null_distribution, perm_diff)
    end
    
    # P-value
    p_perm = mean(abs.(null_distribution) .>= abs(real_diff))
    push!(perm_pvalues_med, p_perm)
end

# FDR correction
adjusted_perm_pvalues_med = adjust(perm_pvalues_med, BenjaminiHochberg())

# Results
top10_idx_med = sortperm(perm_pvalues_med)[1:10]

println("\n=== PERMUTATION TEST RESULTS (MEDIAN) ===")
println("Significant after FDR: $(sum(adjusted_perm_pvalues_med .< 0.05)) / 364")
println("Trend-level (p < 0.005): $(sum(perm_pvalues_med .< 0.005))")

println("\nTop 10 ROIs:")
println("$(lpad("Rank", 5)) $(lpad("ROI Name", 35)) $(lpad("Layer", 15)) $(lpad("Diff", 8)) $(lpad("p-perm", 8)) $(lpad("p-adj", 8))")
println("-"^85)

for (rank, idx) in enumerate(top10_idx_med)
    println("$(lpad(rank, 5)) $(lpad(all_roi_names[idx], 35)) $(lpad(all_roi_layers[idx], 15)) $(lpad(round(perm_differences_med[idx], digits=2), 8)) $(lpad(round(perm_pvalues_med[idx], digits=4), 8)) $(lpad(round(adjusted_perm_pvalues_med[idx], digits=4), 8))")
end

# Create results dataframe
using DataFrames, CSV

results_df = DataFrame(
    ROI_Index = 1:364,
    ROI_Name = all_roi_names,
    Layer = all_roi_layers,
    Median_Diff = perm_differences_med,
    P_Permutation = perm_pvalues_med,
    P_Adjusted_FDR = adjusted_perm_pvalues_med,
    Significant = adjusted_perm_pvalues_med .< 0.05
)

# Sort by p-value
sort!(results_df, :P_Permutation)

# Save
# Save to your analysis folder
output_path = "C:/Users/JLU-MBB/Desktop/Thesis"
CSV.write(joinpath(output_path, "significant_rois_coordinates.csv"), vis_df)
CSV.write(joinpath(output_path, "roi_level_permutation_results.csv"), results_df)

# Create a file with coordinates and statistics for significant ROIs
using DataFrames, CSV

# Get significant + trend ROIs (p < 0.01 for visualization)
visualization_threshold = 0.01
vis_mask = perm_pvalues_med .< visualization_threshold
vis_indices = findall(vis_mask)

println("ROIs to visualize (p < $visualization_threshold): $(length(vis_indices))")

# Separate self and nonself for coordinate extraction
vis_df_list = []

for idx in vis_indices
    if idx <= 37  # Self ROI
        # Parse from self coordinates
        coord_line = self_coords[idx + 1]  # +1 because first line is header
        parts = split(coord_line, "\t")
        
        push!(vis_df_list, (
            ROI_Index = idx,
            ROI_Name = parts[2],
            Layer = parts[6],
            X = parse(Float64, parts[3]),
            Y = parse(Float64, parts[4]),
            Z = parse(Float64, parts[5]),
            Median_Diff = perm_differences_med[idx],
            P_Value = perm_pvalues_med[idx],
            P_Adjusted = adjusted_perm_pvalues_med[idx]
        ))
    else  # Nonself ROI
        # Parse from nonself coordinates
        nonself_idx = idx - 37
        coord_line = nonself_coords[nonself_idx + 1]
        parts = split(coord_line, "\t")
        
        push!(vis_df_list, (
            ROI_Index = idx,
            ROI_Name = parts[2],
            Layer = all_roi_layers[idx],
            X = parse(Float64, parts[3]),
            Y = parse(Float64, parts[4]),
            Z = parse(Float64, parts[5]),
            Median_Diff = perm_differences_med[idx],
            P_Value = perm_pvalues_med[idx],
            P_Adjusted = adjusted_perm_pvalues_med[idx]
        ))
    end
end

vis_df = DataFrame(vis_df_list)

# Save for FSLeyes
CSV.write("significant_rois_coordinates.csv", vis_df)
