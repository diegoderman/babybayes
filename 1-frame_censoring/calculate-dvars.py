from dhcpy import *
import pickle5 as pickle
import multiprocessing as mp
import pandas as pd
import csv
import numpy as np
import argparse


'''
    This script is supposed to run in parallel. The input parameters specific to the dhcp dataset are "hardcoded" both in the dhcpy library and in this script
    The remaining input parameters are as follows:
    input:
 
        --in ~ participants table to use as input, csv or tsv.
    	--out ~ output pickle file path to pickle output.
    	--threads ~ amount of threads to be used.
'''



parser=argparse.ArgumentParser()

parser.add_argument('--threads', type=int, default=48, help='number of threads to use')
parser.add_argument('--in', type=str, default="./results/rel2_full", help='participants table to use as input')
parser.add_argument('--out', type=str, default="./results/rel2_dvars", help='output pickle file path to pickle output-')

# Paths
funcdir = dhcp_root() + "/dhcp_fmri_pipeline"
anatdir = dhcp_root() + "/dhcp_anat_pipeline"

# Participants table to pandas
with open(subjects_pkl + '.pkl', "rb") as fh:
  subjects = pickle.load(fh)
  n_subs = len(subjects.index)

#########################################################################

# get volumes over dvar threshold



core_count = 48

# main parallel function: fsl calculation of dvars.
def f(par, csv, image, bold, boldmask):
    os.system("fsl_motion_outliers -i " + bold + ' -m ' + boldmask + ' -o ' + par + ' --dvars --nomoco -s ' + csv + ' -p ' + image)

# debug parallelization helper
def info(title):
    print(title)
    print('module name:', __name__)
    print('parent process:', os.getppid())
    print('process id:', os.getpid())

# main function
if __name__ == '__main__':
    i = 0
    # initialize incrementers
    thresholds = np.zeros(n_subs, dtype=np.float)
    dvars_means = np.zeros(n_subs, dtype=np.float)
    dvars_sds = np.zeros(n_subs, dtype=np.float)
    
    # run  through the list of subjects 
    while i < len(subjects.index):
        # continue loop if all threads busy
        if (len(mp.active_children())) < core_count - 1:
            
            # Get subject details
            subid = subjects['participant_id'][i]
            sesid = subjects['session_id'][i]

            # Write paths **for dhcp 2nd release, you may need to change this**
            t2mask = anatdir + '/sub-' + subid + "/ses-" + sesid + '/anat/sub-' + subid + '_ses-' + sesid + '_desc-brain_mask.nii.gz'
            boldmask = funcdir + '/sub-' + subid + "/ses-" + sesid + '/func/sub-' + subid + '_ses-' + sesid + '_task-rest_desc-preproc_space-bold_brainmask.nii.gz'
            bold = funcdir + '/sub-' + subid + "/ses-" + sesid + '/func/sub-' + subid + '_ses-' + sesid + '_task-rest_desc-preproc_bold.nii.gz'
            bold2t2w = funcdir + '/sub-' + subid + "/ses-" + sesid + '/xfm/sub-' + subid + '_ses-' + sesid + '_from-bold_to-T2w_mode-image.mat'
            t2w2bold = funcdir + '/sub-' + subid + "/ses-" + sesid + '/xfm/sub-' + subid + '_ses-' + sesid + '_from-T2w_to-bold_mode-image.mat'
            par = funcdir + '/sub-' + subid + "/ses-" + sesid + '/func/sub-' + subid + '_ses-' + sesid + "_desc-DVARS.par"
            csv_fname = funcdir + '/sub-' + subid + "/ses-" + sesid + '/func/sub-' + subid + '_ses-' + sesid + "DVARS.csv"
            image = funcdir + '/sub-' + subid + "/ses-" + sesid + '/func/sub-' + subid + '_ses-' + sesid + "_desc-DVARS-image.png"
            
            # Run main functions
            if not os.path.exists(csv_fname):
                p = mp.Process(target=f, args=(par, csv_fname, image, bold, boldmask))
                p.start()
                print(subid)
            else:
                print("Skipped " + subid + " already has DVARS.")

                # if it already has dvars, accumulate threashold

                # open csv file
                with open(csv_fname) as csv_file:
                    reader = csv.reader(csv_file)
                    individual_dvars = [row[0] for row in reader]

                # make individual dvars into array
                individual_dvars_numpy = np.array(individual_dvars)
                # it defaults to ~float
                float_ind_dvars = individual_dvars_numpy.astype(np.float)
                # get rid of first element always zero
                float_ind_dvars = float_ind_dvars[1:]
                # get threshold 1.5 times iqr over the 75th percentile
                q75, q25 = np.percentile(float_ind_dvars, [75, 25])
                #interquartile range
                iqr = q75 - q25
                #threshold for this individual
                thr_ind = np.percentile(float_ind_dvars, 75) + 1.5*iqr
                # get average dvars
                dvars_ind = np.mean(float_ind_dvars)
                # save in full array of thresholds, means, and std dev
                thresholds[i] = thr_ind
                dvars_means[i] = dvars_ind
                dvars_sds[i] = np.std(float_ind_dvars)
                #np.savetxt('/home/dderman/data.csv', float_ind_dvars, delimiter=',')



        i += 1

    mean_threshold = np.mean(thresholds)
    sd_threshold = np.std(thresholds)
    #print("mean threshold: %f" % mean_threshold)
    print("threshold: %f" % thr_ind)
    print("mean: %f" % dvars_ind)

    print("sd threshold: %f" % sd_threshold)
    print("median dvars: %f" % np.median(dvars_means))
    print("mean sd: %f" % np.mean(dvars_sds))
