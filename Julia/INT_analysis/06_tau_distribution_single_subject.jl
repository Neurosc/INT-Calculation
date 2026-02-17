
using Plots

# Some colors and fonts to use later on 
pltcolors= cgrad(:matter, 5, categorical = true)
default(fontfamily="Computer Modern")  # LaTeX-style font

# Load data for one subject
subject_id = "048"
self_file = "C:/Users/JLU-MBB/Desktop/Thesis/ACW_results/controls/timeseries_extracted_self_sub-$subject_id.jld2"
nonself_file = "C:/Users/JLU-MBB/Desktop/Thesis/ACW_results/controls/timeseries_extracted_nonself_sub-$subject_id.jld2"



# Load both files
self_data = load(self_file)
nonself_data = load(nonself_file)

# Get first scan from each
self_scan = self_data["scans"][1]
nonself_scan = nonself_data["scans"][1]

# Get tau values to find representative ROIs
self_tau = self_scan.acw_results[1]
nonself_tau = nonself_scan.acw_results[1]

# Select 3 representative self ROIs (low, medium, high tau)
self_tau_sorted_idx = sortperm(self_tau)
self_roi_low = self_tau_sorted_idx[round(Int, length(self_tau) * 0.25)]
self_roi_med = self_tau_sorted_idx[round(Int, length(self_tau) * 0.50)]
self_roi_high = self_tau_sorted_idx[round(Int, length(self_tau) * 0.75)]

# Select 3 representative nonself ROIs (low, medium, high tau)
nonself_tau_sorted_idx = sortperm(nonself_tau)
nonself_roi_low = nonself_tau_sorted_idx[round(Int, length(nonself_tau) * 0.25)]
nonself_roi_med = nonself_tau_sorted_idx[round(Int, length(nonself_tau) * 0.50)]
nonself_roi_high = nonself_tau_sorted_idx[round(Int, length(nonself_tau) * 0.75)]

# Get lags
lags = self_scan.lags

# Create plot
p = plot(xlabel="Lag (seconds)",
         ylabel="Autocorrelation",
         title="Autocorrelation Function of Sub-$subject_id",
         legend=:topright,
         size=(900, 600),
         legendfontsize=9)

# Plot self ROIs (blue shades)
plot!(lags, self_scan.acf[self_roi_low, :],
      label="Self ROI $self_roi_low (τ=$(round(self_tau[self_roi_low], digits=2))s)",
      color=pltcolors[2], linewidth=2, linestyle=:solid)

plot!(lags, self_scan.acf[self_roi_med, :],
      label="Self ROI $self_roi_med (τ=$(round(self_tau[self_roi_med], digits=2))s)",
      color=pltcolors[2], linewidth=2, linestyle=:solid)

plot!(lags, self_scan.acf[self_roi_high, :],
      label="Self ROI $self_roi_high (τ=$(round(self_tau[self_roi_high], digits=2))s)",
      color=pltcolors[2], linewidth=2, linestyle=:solid)

# Plot nonself ROIs (red shades)
plot!(lags, nonself_scan.acf[nonself_roi_low, :],
      label="Nonself ROI $nonself_roi_low (τ=$(round(nonself_tau[nonself_roi_low], digits=2))s)",
      color=pltcolors[4], linewidth=2, linestyle=:solid)

plot!(lags, nonself_scan.acf[nonself_roi_med, :],
      label="Nonself ROI $nonself_roi_med (τ=$(round(nonself_tau[nonself_roi_med], digits=2))s)",
      color=pltcolors[4], linewidth=2, linestyle=:solid)

plot!(lags, nonself_scan.acf[nonself_roi_high, :],
      label="Nonself ROI $nonself_roi_high (τ=$(round(nonself_tau[nonself_roi_high], digits=2))s)",
      color=pltcolors[4], linewidth=2, linestyle=:solid)

# Add zero line
hline!([0], linestyle=:dot, color=:black, linewidth=1, label="")

display(p)


