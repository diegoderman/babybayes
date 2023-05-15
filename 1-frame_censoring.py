#!/usr/bin/env python
#
# Diego Derman
# FRG - IUB


"""
    filter-dvars.py: Generates dHCP participants list after frame censoring from DVARS threshold.
    
    Description: This scripts generates a new output table with the amount of frames over the subject-specific DVARS threshold (1.5 times over IQR). It is also used to filter out the subjects with too many frames over the subject-specific DVARS threshold (less than 1600 remianing frames).
    The main function takes the previously created DVARS table in the canonical DHCP func subject directory, (it must have been previously created using the `fsl_motion_outliers` command).
    The script is written to run in a multithread setting. The input parameters specific to the dhcp dataset are "hardcoded" both in the dhcpy library and in this script
    
    Input: The remaining input parameters are as follows:
 
        --in ~ participants table to use as input, csv or tsv.
    	--out ~ output pickle file path to pickle output. After DVARS frame censoring
    	--threads ~ amount of threads to be used.
"""

########################################################################
# Dependencies

from dhcpy import *
import pickle5 as pickle
import multiprocessing as mp
import pandas as pd
import csv
import numpy as np
import argparse

#########################################################################
# Arguments and hard-coded (dataset-specific) parameteres

# CLI Argument parser 

parser=argparse.ArgumentParser()
parser.add_argument('--threads', type=int, default=48, help='number of threads to use')
parser.add_argument('--in', type=str, default="./results/rel2_full", help='participants table to use as input')
parser.add_argument('--out', type=str, default="./results/rel2_dvars", help='output pickle file path to pickle output-')

# Assign parsed variables
args=parser.parse_args()

core_count = args.threads
subjects_pkl = args.in
subjects_dvars = args.out

# DHCP Paths
funcdir = dhcp_root() + "/dhcp_fmri_pipeline" # dhcp_root() is a function of dhcpy
anatdir = dhcp_root() + "/dhcp_anat_pipeline"

#########################################################################
# Parallelization functions and main thread


# Participants table to pandas
with open(subjects_pkl + '.pkl', "rb") as fh:
  subjects = pickle.load(fh)


# main function
def dvars(csvfile):

    '''
	 dvars: calculates the interval of minimum amount of frames over the DVARS threshold as explained in Eyre et al., 2022. Returns the amount of frames over the threshold, the index where the interval starts, and the mean dvars over the whole run. 
	 
	 Input: The previously generated csvfile of a single subject's dvars as input (see fsl_motion_outliers).
	 
         returns the rating and location of subjects 1600 best volumes
         returns -1 if the csvfile cannot be found.
    '''
    
    # if file does not exist, return -1
    if not os.path.exists(csvfile):
        print(csvfile)
        return -1
    else:
        # else, open csv file as pandas and format correctly
        subject_dvars = csv2pd(csvfile)
        subject_dvars = subject_dvars.astype(float)
        # save mean dvars over all the frames
        mean_dvars = np.mean(subject_dvars[1:])[0]
        
        # definition of motion outlier **for dhcp 2nd release, you may need to change this**
        min_volsOver = 1600
        keep_vols = 1600
        BOLD_length = 2299
       
        # definition of threshold as > 1.5 IQR over 75th percentile of each subject's DVARS. 
        threshold = 1.5 * stats.iqr(subject_dvars[1:]) + np.percentile(subject_dvars[1:], 75)

        # Intialize index of lowest motion interval
        min_idx = 0
        for i in range (BOLD_length-keep_vols):
            # Intialize value of lowest motion interval at 1600
            previous_min = min_volsOver
            # if the moving window contains less than 1600 motion outlier frames, keep it
            min_volsOver = min(min_volsOver, np.count_nonzero(subject_dvars[i:i+keep_vols] > threshold))
            # if the value was not kept, also save the new index of the window.
            if previous_min != min_volsOver:
                    min_idx = i
    return [min_volsOver, min_idx, mean_dvars]

# main thread

# initialize rows of final table
rows = list()

for participant in subjects.iterrows():
    # get first list of the row (format)
    participant=participant[1]
    # get participant information for paths
    subid = participant['participant_id']
    sesid = str(participant['session_id'])
    # dhcp canonical paths
    session_path = funcdir + "/sub-" + subid + "/ses-" + sesid
    fd_path = session_path + '/func/sub-' + subid + '_ses-' + sesid + '_motion.tsv'
    subject_fd = csv2pd(fd_path)
    subject_fd = subject_fd.astype(float)
    # calculate mean framewise displacement for the subject.
    mean_fd = subject_fd['framewise_displacement'].mean()
    # if the DVARS file was previously created, then calculate DVARS outlier interval.
    if os.path.exists(session_path + '/func/sub-' + subid + '_ses-' + sesid + 'DVARS.csv'):
        # dvars main function
        vols_dvars, index_dvars, mean_dvars = dvars(session_path + '/func/sub-' + subid + '_ses-' + sesid + 'DVARS.csv')
        
        # append subject data and output of dvars to the final table.
        rows.append({'participant_id': subid,
                     'singleton': participant['singleton'],
                     'birth_age': float(participant['birth_age']),
                     'sex': participant['sex'],
                     'birth_weight': float(participant['birth_weight']),
                     'session_id': participant['session_id'],
                     'scan_age': float(participant['scan_age']),
                     'scan_number': int(participant['scan_number']),
                     'dvars_outliers': vols_dvars,
                     'mean_dvars': mean_dvars,
                     'start_best_interval': index_dvars,
                     'mean_fd': mean_fd
                     })
    else:
        print("warning: dvars does not exist. sub: " + subid)
        
# discard first (empty) row.
rows.pop(0)


# tidy final table
subjects_out = pd.DataFrame(rows, columns=['participant_id', 'singleton', 'birth_age', 'sex', 'birth_weight',
                                           'session_id', 'scan_age', 'scan_number', 'dvars_outliers', 'mean_dvars',
                                           'start_best_interval', 'mean_fd'])

# Save in pkl and csv on subject_dvars
subjects_out.to_pickle(subjects_dvars + '.pkl')
subjects_out.to_csv(subjects_dvars + '.csv')
