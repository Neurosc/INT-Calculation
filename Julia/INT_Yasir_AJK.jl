#This code is for calculating the intrinsic Neural timescales of 
#Glaser-parcellated data.

using IntrinsicTimescales, MAT, Statistics, Plots

#parameters
TR = 2.0 #seconds
fs = 1/TR #sampling interval

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
                acwtypes=[:auc], # tau and tau will work the best for my data probably...
                skip_zero_lag = false, # fs is too low, already goes down significantly after first lag so will only make noise
                n_lags = 40) # Not sure if this many lags are too much... I have 90 tps
            push!(all_trials_data, result)
        end 
    end
    

    return all_trials_data

end

base_path = "C:/Users/JLU-MBB/Desktop/MYELIN-main/DATA/qin_organized"
groups = ["controls", "meditators"]

# Initialize dictionaries to store results by subject
controls_subjects = Dict()
meditators_subjects = Dict()

for group in groups
    group_path = joinpath(base_path, group)
    subjects = readdir(group_path, join=false)
    subjects = filter(s -> startswith(s, "sub-"), subjects)

    for subject in subjects 
        subject_path = joinpath(group_path, subject)
        all_scans = process_subject(subject_path)  # Returns array of scan results for this subject
        
        # Store all scans for this subject together
        if group == "controls"
            controls_subjects[subject] = all_scans  # Key: subject name, Value: array of all their scans
        else
            meditators_subjects[subject] = all_scans
        end
    end
end

using JLD2

# Create output directories
output_dir = "C:/Users/JLU-MBB/Desktop/MYELIN-main/DATA/ACW_results"
controls_dir = joinpath(output_dir, "controls")
meditators_dir = joinpath(output_dir, "meditators")

mkpath(controls_dir)
mkpath(meditators_dir)

# Save each subject as a separate file
for (subject_id, scans) in controls_subjects
    @save joinpath(controls_dir, "$(subject_id).jld2") scans
end

for (subject_id, scans) in meditators_subjects
    @save joinpath(meditators_dir, "$(subject_id).jld2") scans
end

println("Saved $(length(controls_subjects)) controls subjects")
println("Saved $(length(meditators_subjects)) meditators subjects")