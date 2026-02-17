# Load one file
test_file = "C:/Users/JLU-MBB/Desktop/ACW_results/controls/timeseries_extracted_self_sub-048.jld2"
data = load(test_file)
scan1 = data["scans"][1]

println("acwtypes: ", scan1.acwtypes)
println("acw_results type: ", typeof(scan1.acw_results))
println("acw_results size/length: ", size(scan1.acw_results))
println("Tau values (first 5): ", scan1.acw_results[1][1:5])
println("AUC values (first 5): ", scan1.acw_results[2][1:5])
println("Number of ROIs in tau: ", length(scan1.acw_results[1]))
println("Number of ROIs in auc: ", length(scan1.acw_results[2]))


using Plots

# Get ACF for one ROI from one scan
scan1 = data["scans"][1]

# ACF is stored in scan1.acf (ROIs × lags)
roi_1_acf = scan1.acf[5, :]  # First ROI's autocorrelation across 40 lags

# Get the lags (in seconds)
lags = scan1.lags  # Should be 0, 2, 4, 6, ... (TR = 2 seconds)

# Plot
plot(lags, roi_1_acf, 
     xlabel="Lag (seconds)",
     ylabel="Autocorrelation",
     title="ACF for ROI 5",
     label="ACF",
     linewidth=2,
     marker=:circle)

# Add horizontal line at zero
hline!([0], linestyle=:dash, color=:black, label="Zero")

println("Tau for ROI 1: ", scan1.acw_results[1][1])
println("AUC for ROI 1: ", scan1.acw_results[2][1])


# Find an ROI with longer timescale
tau_values = scan1.acw_results[1]
longest_roi = argmax(tau_values)

println("\nROI with longest tau (ROI $longest_roi):")
println("  Tau: ", tau_values[longest_roi], " seconds")

# Plot it
plot(lags, scan1.acf[longest_roi, :],
     xlabel="Lag (seconds)",
     ylabel="Autocorrelation", 
     title="ACF for ROI $longest_roi (longest tau)",
     label="ACF",
     linewidth=2,
     marker=:circle)
hline!([0], linestyle=:dash, color=:black, label="Zero")

println("ROI 5:")
println("  Tau: ", scan1.acw_results[1][5], " seconds")
println("  AUC: ", scan1.acw_results[2][5])

println("\nROI 7 (longest):")
println("  Tau: ", scan1.acw_results[1][7], " seconds")  
println("  AUC: ", scan1.acw_results[2][7])
