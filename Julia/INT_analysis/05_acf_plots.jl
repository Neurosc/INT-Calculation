# This file visualizes autocorrelation functions for representative self and nonself ROIs. 


# Plots illustrate:
#### the range of intrinsic timescales across brain regions for a single subject.



# Determine common axis limits
max_tau = max(maximum(all_self_tau), maximum(all_nonself_tau))

# Create overlapping histogram
p = histogram(all_self_tau,
              bins=50,
              alpha=0.6,
              label="Self ROIs (N=$(length(all_self_tau)))",
              color=:blue,
              xlims=(0, max_tau),
              xlabel="Tau (seconds)",
              ylabel="Count",
              title="Distribution of Tau Values Across All Subjects",
              size=(900, 600))

histogram!(all_nonself_tau,
           bins=50,
           alpha=0.6,
           label="Nonself ROIs (N=$(length(all_nonself_tau)))",
           color=:red)

# Add mean lines
vline!([mean(all_self_tau)], linewidth=2, color=:blue, linestyle=:dash, label="Self Mean")
vline!([mean(all_nonself_tau)], linewidth=2, color=:red, linestyle=:dash, label="Nonself Mean")

display(p)
savefig("tau_distribution_all_subjects.png")



using Plots, Statistics, StatsBase


# CREATE TWO SEPARATE HISTOGRAMS WITH PERCENTAGE


# Determine common axis limits
max_tau = max(maximum(all_self_tau), maximum(all_nonself_tau))

# Create histograms and manually calculate percentages
edges = 0:2.5:max_tau  # Bin edges

hist_self = fit(Histogram, all_self_tau, edges)
hist_nonself = fit(Histogram, all_nonself_tau, edges)

# Convert counts to percentages
pct_self = (hist_self.weights ./ length(all_self_tau)) .* 100
pct_nonself = (hist_nonself.weights ./ length(all_nonself_tau)) .* 100

max_pct = max(maximum(pct_self), maximum(pct_nonself))

# Self ROIs histogram
p1 = bar(hist_self.edges[1][1:end-1],
         pct_self,
         title="Self ROIs (N=$(length(all_self_tau)))",
         xlabel="Tau (seconds)",
         ylabel="Percentage (%)",
         label="",
         alpha=0.7,
         color=pltcolors[3],
         xlims=(0, max_tau),
         ylims=(0, max_pct * 1.1),
         bar_width=2.5)
vline!([mean(all_self_tau)], linewidth=2, color=:darkblue, linestyle=:dash, label="Mean")
vline!([median(all_self_tau)], linewidth=2, color=:orange, linestyle=:dot, label="Median")

# Nonself ROIs histogram
p2 = bar(hist_nonself.edges[1][1:end-1],
         pct_nonself,
         title="Nonself ROIs (N=$(length(all_nonself_tau)))",
         xlabel="Tau (seconds)",
         ylabel="Percentage (%)",
         label="",
         alpha=0.7,
         color=pltcolors[4],
         xlims=(0, max_tau),
         ylims=(0, max_pct * 1.1),
         bar_width=2.5)
vline!([mean(all_nonself_tau)], linewidth=2, color=:darkred, linestyle=:dash, label="Mean")
vline!([median(all_nonself_tau)], linewidth=2, color=:orange, linestyle=:dot, label="Median")

# Combine into single figure
p = plot(p1, p2,
         layout=(2,1),
         size=(900, 800),
         plot_title="Distribution of Tau Values Across All Subjects")

display(p)
savefig("tau_distribution_percentage_all_subjects.png")