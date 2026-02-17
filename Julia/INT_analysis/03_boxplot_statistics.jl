# This file performs statistical comparisons between self and nonself ROIs within and between groups. 
#### Results are visualized as boxplots.

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


# PRINT STATISTICS SUMMARY

println("\n=== STATISTICAL SUMMARY ===")
println("\nWithin-Group (Wilcoxon signed-rank):")
println("  Controls (Self vs Nonself): p = ", round(p_control_paired, digits=4))
println("  Meditators (Self vs Nonself): p = ", round(p_med_paired, digits=4))

println("\nBetween-Group (Mann-Whitney U):")
println("  Self ROIs (C vs M): p = ", round(p_self_between, digits=4))
println("  Nonself ROIs (C vs M): p = ", round(p_nonself_between, digits=4))

println("\n* p < 0.05, ** p < 0.01, *** p < 0.001")


# Boxplot for meditators vs controls

control_labels = fill("Controls", length(control_all))
meditator_labels = fill("Meditators", length(meditator_all))

all_data_groups = vcat(control_all, meditator_all)
all_labels_groups = vcat(control_labels, meditator_labels)

# Get p-value for annotation
p_overall = pvalue(overall_test)
y_max = maximum(all_data_groups)

# Create boxplot
p = boxplot(all_labels_groups, all_data_groups,
        ylabel="Intrinsic Timescale (seconds)",
        title="Overall Meditation Effect on Intrinsic Timescales",
        legend=false,
        fillalpha=0.7,
        linewidth=2,
        size=(600, 700),
        color=pltcolors[3],
        ylims=(0, y_max * 1.15))

# Add individual points
dotplot!(all_labels_groups, all_data_groups,
         marker=(:circle, 4, 0.4),
         color=pltcolors[1])

# Add significance bar
plot!([1, 2], [y_max * 1.05, y_max * 1.05], color=:black, linewidth=1.5)
annotate!(1.5, y_max * 1.08,
         text(p_overall < 0.001 ? "***" : p_overall < 0.01 ? "**" : p_overall < 0.05 ? "*" : "ns",
              10, :center))

display(p)
savefig("overall_meditation_effect.png")