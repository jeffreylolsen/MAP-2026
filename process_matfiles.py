import sys
import site

site.addsitedir(site.getusersitepackages())
import scipy
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sklearn
from scipy.io import loadmat
import sktime
import pyarrow

## Directory of all .MAT files
matdata_dir = "./data/matfiles"

## File with variable details
varfile = pd.read_csv("./matfile_variables.csv")

## Disposition
disposition = pd.read_csv("./data/disposition.csv")

## Subject drive info (including fatigue)
subj_drives = pd.read_csv("./data/subject_drive_info.csv")

## All analyses should only use files marked "Reduced"
reduced_disp = disposition[disposition["Reduced"] == "X"]

## Add .mat to DAQs so they can be read from the mat data directory, reset index to facilitate looping
reduced_disp = reduced_disp.assign(
    MatName=reduced_disp["DaqName"].str[:-4] + ".mat"
).reset_index(0, drop=True)

## Ensure order of columns
final_df = pd.DataFrame(columns=varfile["Ordered Columns"].tolist())

# Read in each mat file and
for i in range(0, reduced_disp.shape[0] - 1):

    ## Read current mat file
    fname = matdata_dir + "/" + reduced_disp["MatName"][i]
    try:
        data = loadmat(fname)
        print(f"Successfully read {fname} ({i + 1}/{reduced_disp.shape[0] - 1})")
    except FileNotFoundError:
        print(f"Skipping {fname}: File not found ({i + 1}/{reduced_disp.shape[0] - 1})")

    # Names of variables to extract, order does not matter since final_df is initialized
    variable_names = dict(
        zip(varfile["Variable Names"].dropna(), varfile["Original Names"].dropna())
    )
    variable_series = {}

    ## Full lat pos (of vehicle) series (shaped differently so taken differently)
    ## And road width to compare to lane deviation
    variable_series["Vehicle_Lat_Dev"] = (
        data["elemDataI"]["SCC_Lane_Deviation"][0][0][:, 1]
        .flatten()
    )
    variable_series["Road_Width"] = (
        data["elemDataI"]["SCC_Lane_Deviation"][0][0][:, 2]
        .flatten()
    )
    
    # Code task availability from LogStream
    variable_series["Task_Available_raw"] = (
        data["elemDataI"]["SCC_LogStreams"][0][0][:, 2]
        .flatten().astype(int)
    )
    variable_series["Task_Available"] = variable_series["Task_Available_raw"].astype(str)
    variable_series["Task_Available"][variable_series["Task_Available_raw"] <= 1] = "None"
    variable_series["Task_Available"][variable_series["Task_Available_raw"] == 2] = "Distractor"
    variable_series["Task_Available"][variable_series["Task_Available_raw"] == 4] = "Left"
    variable_series["Task_Available"][variable_series["Task_Available_raw"] == 5] = "Right"
    variable_series.pop("Task_Available_raw")

    # Get number of frames
    nstep = variable_series["Vehicle_Lat_Dev"].shape[0]

    for key, val in variable_names.items():
        try:
            variable_series[key] = data["elemDataI"][val][0][0].flatten()
        except ValueError:
            variable_series[key] = np.repeat(99, nstep)

    ## Assemble DF w/ Subject and Drive meta data
    temp_df = pd.DataFrame(variable_series)
    temp_df.insert(0, "Subject", np.repeat(reduced_disp["Subject"][i], nstep))
    temp_df.insert(1, "Drive", np.repeat(reduced_disp["DriveN"][i], nstep))
    temp_df["DaqName"] = np.repeat(reduced_disp["MatName"][i][:-4], nstep)

    ## Attach "Sample" for each 60 sec sample
    temp_df["Sample"] = 1 + temp_df.index // 3600

    ## Discard any samples that aren't length 3600
    # NOTE might not want to discard these if we aren't using samples
    # NOTE samples are not all equal lengths if below line is commented out
    # temp_df = temp_df[temp_df.groupby("Sample")["Sample"].transform("size") == 3600]

    ## Append the temp df to the current final df
    final_df = pd.concat([final_df, temp_df], ignore_index=True)

# Convert eyelid closed variables to boolean
print("Processing eyelid closed variables...")
for side in ["Left", "Right"]:
    final_df[f"{side}_Eyelid_Closed"] = (
        final_df[f"{side}_Eyelid_Closed"]
        .apply(lambda x: str(x[0]).lower() if isinstance(x, list) else str(x).lower())
        .replace({"['true']": True, "['false']": False, "99": np.nan, "[]": np.nan})
        .astype("boolean")
    )

# Add fatigue and BAC info
# NOTE NA values in some cases on restarts because they weren't asked again before the second Daq, currently this is unhandled
print("Processing fatigue, RTD, and BAC info...")
subj_drives["DaqName"] = subj_drives["Daq"].str[:-4]
subj_drives = subj_drives.loc[
    :, ["DaqName", "rtd", "csra_kss_score", "start_BAC"]
].rename(
    columns={
        "rtd": "Ready_to_Drive",
        "csra_kss_score": "KSS_Score",
        "start_BAC": "Start_BAC",
    }
)
subj_drives["Ready_to_Drive"] = subj_drives["Ready_to_Drive"] == 1
final_df = pd.merge(final_df, subj_drives, on="DaqName", how="left")

# Generate sample IDs
print("Processing sample IDs...")
final_df["Sample_ID"] = "ID_" + final_df.loc[:, ["Subject", "Drive", "Sample"]].astype(
    str
).agg("_".join, axis=1)

# Make factor variables
print("Converting to factor variables...")
final_df["DaqName"] = final_df["DaqName"].astype("category")
final_df["Subject"] = final_df["Subject"].astype("category")
final_df["Drive"] = final_df["Drive"].astype("category")
final_df["Sample_ID"] = final_df["Sample_ID"].astype("category")
final_df["Task_Available"] = final_df["Task_Available"].astype("category")

print("Writing output...")
final_df.to_feather("./data/non_aug_data.feather")
