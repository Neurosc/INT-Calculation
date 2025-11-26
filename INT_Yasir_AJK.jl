#This code is for calculating the intrinsic Neural timescales of 
#Glaser-parcellated data.


using IntrinsicTimescales, MAT, Statistics, Plots

#parameters
TR = 2.0 #seconds
fs = 1/TR #sampling interval

# Initialize storage arrays for each group
# Initialize storage arrays for each group
control_auc = []
control_tau = []
control_results = []
meditator_auc = []
meditator_tau = []
meditator_results = []

#Load region classification
mat_data = matread("C:/Users/JLU-MBB/Desktop/MYELIN-main/DATA/myDataParcels.mat")
parcel_class = vec(mat_data["myDataParcels"][:,2]) # from mydataparcels, get the second column
self_mask = parcel_class .== 2 # 33 self regions
nonself_mask = parcel_class .== 1 # 327 nonself regions

# Function to process one subject
# Each rep and each meditation type is processed separately.
function process_subject(subject_path)

    scan_types = ["ca1", "ca2", "conc1", "conc2", "metta1", "metta2"]

    #Storage for this subject's data
    all_trials_data = []

    for scan in scan_types
        
        filename = joinpath(subject_path, "$(basename(subject_path))_$(scan).mat")
        
        if isfile(filename)
            mat = matread(filename)

            timeseries = mat["parcellated"] # this gives parcels x timepoints
            result= acw(timeseries, fs,
                dims=2, #the timeseries is in 2nd dim
                acwtypes=[:auc, :tau], # tau and tau will work the best for my data probably...
                skip_zero_lag = true, # fs is too low, already goes down significantly after first lag so will only make noise
                n_lags = 40) # Not sure if this many lags are too much... I have 90 tps
            push!(all_trials_data, result)
        end 
    end
    

    return all_trials_data

end

base_path = "C:/Users/JLU-MBB/Desktop/parcellated"
groups = ["control", "meditator"]

for group in groups # Loop for the categories
    group_path = joinpath(base_path, group)
    subjects = readdir(group_path, join=false)
    subjects = filter(s -> startswith(s, "sub-"), subjects)

    for subject in subjects 
        subject_path = joinpath(group_path, subject)
        all_scans = process_subject(subject_path) # Feed it in subject by subject at each category

        # Store each scan's ACW results
        for scan_result in all_scans
            tau = scan_result.acw_results[2]  # 360 values
            auc = scan_result.acw_results[1]

            if group == "control"
                push!(control_auc, auc)
                push!(control_tau, tau)
                push!(control_results, scan_result)  # Store full object

            else
                push!(meditator_auc, auc)
                push!(meditator_tau, tau)
                push!(meditator_results, scan_result)  # Store full object

            end
        end
    end
end

#############################################
#Some Test Plots...
#############################################


#### 1- Give the ACW curve of the first subject ####

#  one scan to validate
test_scan = control_results[1]

# Get the data
lags = collect(test_scan.lags)  # Convert range to array
acf_data = test_scan.acf  # 360 parcels × 30 lags
auc_values = test_scan.acw_results[1]  # AUC values
tau_values = test_scan.acw_results[2]  # TAU values (fitted τ)


# Choose first 10 parcels to plot
parcels_to_plot = 1:10 

# Create plot
p = plot(title="ACW of Some Parcels (tau)",
         xlabel="Lag (seconds)",
         ylabel="Autocorrelation",
         legend=:topright,
         size=(800, 600))

# Plot ACF curve for each parcel
for parcel in parcels_to_plot
    acf_val = acf_data[parcel, :]  # Get ACF for this parcel across all lags
    plot!(p, lags, acf_val, 
          label="Parcel $parcel",
          linewidth=2,
          alpha=0.7)
end

# Add zero line
hline!(p, [0], color=:black, linestyle=:dash, linewidth=2, label="Zero")
display(p)

#Print the values given for that acw function
tau_values = test_scan.acw_results[2] 
println("TAU values for 10 parcels:")
for i in 1:20
    println("  Parcel $i: tau = $(round(tau_values[i], digits=3))") # , TAU = $(round(tau_values[i], digits=3))")
end


#### T-Test for Self vs Nonself (Control) ####

using HypothesisTests, Statistics

println("\n" * "="^60)
println("CONTROL GROUP: Self vs Non-Self")
println("="^60)

# Convert to matrix
control_tau_matrix = hcat(control_tau...)'  # trials × 360 parcels

# Pool all trials together - flatten into one long vector
control_self_all = vec(control_tau_matrix[:, self_mask])  # All self values
control_nonself_all = vec(control_tau_matrix[:, nonself_mask])  # All non-self values

# Perform t-test
test_result = UnequalVarianceTTest(control_self_all, control_nonself_all)


#### 2- Fit exponential decay and compare with empirical ACF(tau) ####
plots_array = []

parcels = [1,2,3,4,5,6,7,8,9,10]
for parcel in parcels 
    #Get the empirical acf
    empirical_acf = acf_data[parcel,:]

    #Get the fitted tau
    tau_fitted = tau_values[parcel]

    # generate an exponential function with the fitted tau
    theoretical_acf = exp.(-lags ./ tau_fitted)

    # Region type control
    region_type = "control"

    # Create plot
    p_comp = plot(title="Parcel $parcel ($region_type)",
                  xlabel="Lag (seconds)",
                  ylabel="Autocorrelation",
                  legend=:topright)
    
                      # Plot empirical ACF
    plot!(p_comp, lags, empirical_acf,
          label="Empirical ACF",
          linewidth=3,
          color=:blue,
          alpha=0.7)
    # Plot theoretical exponential
    plot!(p_comp, lags, theoretical_acf,
          label="exp(-l/τ)",
          linewidth=2,
          linestyle=:dash,
          color=:red)
        # Zero line
    hline!(p_comp, [0], color=:black, linestyle=:dot, linewidth=1, label="")

        # Calculate goodness of fit (R²)
    ss_res = sum((empirical_acf - theoretical_acf).^2)
    ss_tot = sum((empirical_acf .- mean(empirical_acf)).^2)
    r_squared = 1 - ss_res / ss_tot
    
    # Annotate
    annotate!(p_comp, lags[end]*0.6, 0.9,
             text("τ = $(round(tau_fitted, digits=2))s\nR² = $(round(r_squared, digits=3))",
                  :left, 8))
    
    push!(plots_array, p_comp)
    
    println("Parcel $parcel ($region_type): τ = $(round(tau_fitted, digits=3))s, R² = $(round(r_squared, digits=3))")
end

# Display grid
plot(plots_array..., layout=(2, 3), size=(1400, 900))



