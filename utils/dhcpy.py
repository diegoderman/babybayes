import os
import csv
import pandas as pd
import numpy as np
from scipy import stats
from socket import gethostname
import subprocess
from socket import gethostname
import pickle5 as pickle

'''
Python tools for processing the dHCP dataset. Most functions 
would also be useful for any other fsl pipeline.
'''

def get_host():
    hostname = gethostname()
    if hostname != "viper":
        if hostname[0:2] == "dl" or hostname[0] == 'i':
            return "job_manager"
        else: return "carbonate"
    else: return "viper"

def init_host():
    carbonate_home = "/N/u/dderman/Carbonate"
    baby_ICA = "/N/project/baby_ICA"
    core_count = 24
    host = get_host()
    if host == "job_manager":
        carbonate_home = "/N/u/dderman/Carbonate"
        baby_ICA = "/N/project/baby_ICA"
        core_count = 24
    elif host == "carbonate":
        carbonate_home = "/N/u/dderman/Carbonate"
        baby_ICA = "/N/project/baby_ICA"
        core_count = 46
    elif host == "viper":
        carbonate_home = "/home/dderman/carbonate"
        baby_ICA = "/home/dderman/baby_ICA"
        core_count = 7
    else:
        print("Warning: unrecognized host, trying carbonate paths.")
        [carbonate_home, baby_ICA, core_count] = [-1, -1, -1]

    return [carbonate_home, baby_ICA, core_count]

def dhcp_root():
    '''
     string: Returns the full path to the root directory for the dhcp dataset.
    '''
    if gethostname() == 'viper':
        return "/home/dderman/baby_ICA"
    else:
	    return "/N/project/baby_ICA"

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def success_report(dataframe, output_field, output, subid, sesid):
    '''
        To the dataframe, adds a column stating the output state of the algorithm.
    '''

    #result = dataframe

    if output_field not in dataframe.columns:
        return 1
    else:
        dataframe.loc[dataframe['session_id'] == sesid, 'apply_mask'] = output
        #result.loc[result['session_id'] == sesid, 'apply_mask'] = output
        return 0



def is_multisession(dataframe, subid):
    '''
        input: pandas dataframe. DEPRECATED FROM LAST VERSION NO LONGER A PKLFILE
        returns if a subject has multiple session in said dataframe.
    '''

    sessions = dataframe.loc[dataframe['participant_id'] == subid]
    amount_ses = len(sessions)
    if amount_ses == 1:
        return 0
    else:
        return sessions['session_id'].tolist()

def filter_by(input_dataframe, expression) -> object:
    '''
        input: pandas dataframe. DEPRECATED FROM LAST VERSION NO LONGER A PKLFILE
        returns a filtered pandas of the csvfile list with the n lowest (highest) after sorted by field.
        expression: must be a statement of the form field >,<,= value.
    '''

    # parse to check if spaces were added
    import re

    if " " not in expression:
        '''
        split = re.split('>|<|=|>=|<=|==', expression)
        field = split[0]
        value = split[1]
        symbol = expression[len(field) + 1:len(field) + 2]
        '''
        print("ERROR: Expression must contain spaces.")
        return 1
    else:
        split = expression.split()
        field = split[0]
        symbol = split[1]
        value = split[2]
        if is_number(value): value = float(value)

    # load pandas dataframe
    df = input_dataframe

    # Make sure field is in the header

    if field not in df.columns:
        print("ERROR: Unknown field " + field + ".")
        print(df.columns)

    # get type of field.

    if type(df[field]) == str:
        print("WARNING: Field " + field + " was str.")
        df[field] = df[field].to_numeric()

    if symbol == '>':
        return df[df[field] > value]
    elif symbol == '<':
        return df[df[field] < value]
    if symbol == '>=':
        return df[df[field] >= value]
    if symbol == '<=':
        return df[df[field] <= value]
    if symbol == '=':
        return df[df[field] == value]
    if symbol == '==':
        return df[df[field] == value]
    else:
        print("Error symbol.")
        return 1

def dvars(csvfile):
    '''
		 csvfile of a single subject's dvars as input
         returns the rating and location of subjects 1600 best volumes
         returns -1 is the subject is to be discarded according to dvars.
    '''
    if not os.path.exists(csvfile):
        print(csvfile)
        return -1
    else:
        subject_dvars = csv2pd(csvfile)
        subject_dvars = subject_dvars.astype(float)
        mean_dvars = np.mean(subject_dvars[1:])[0]
        #zscore = stats.zscore(subject_dvars['0']) wrong
        min_volsOver = 1600
        #threshold = 1.5*stats.iqr(zscore) + np.percentile(zscore, 75) wrong
        threshold = 1.5 * stats.iqr(subject_dvars[1:]) + np.percentile(subject_dvars[1:], 75)

        min_idx = 0
        for i in range (2299-1600):
            premin = min_volsOver
            #print(np.count_nonzero((subject_dvars[i:i+1600]) > threshold))
            min_volsOver = min(min_volsOver, np.count_nonzero(subject_dvars[i:i+1600] > threshold))
            if premin != min_volsOver:
                    min_idx = i
    return [min_volsOver, min_idx, mean_dvars]

def get_niftiPath(participant_id, session="bold"):
    '''
    :param participant_id: a string containing the participant id
    :param session: a string either 'bold' or 'dmri', 'anat' NOT IMPLEMENTED
    :return niftiPath: a string containing the full boldPath
    returns false if no sessions or multiple sessions.
    '''

    dhcp_dir = '/geode2/home/u020/dderman/Carbonate/dHCP_dataset/'
    #dhcp_dir = '/N/project/baby_ICA/'
    nosessions_flag = False
    nopart_flag = False

    # Implement switch case
    if session == "bold":
        subject_path = dhcp_dir + 'dhcp_fmri_pipeline/sub-' + participant_id
    else:
        subject_path = dhcp_dir + 'dhcp_fmri_pipeline/sub-' + participant_id
    amount_sessions = 0
    if os.path.exists(subject_path):
        for directory in os.listdir(subject_path):
            if directory[0:3] == 'ses':
                amount_sessions = amount_sessions + 1
                session_name = directory
        if amount_sessions < 1:
            nosessions_flag = True
            niftiPath = False
        else:
            niftiPath = subject_path + '/' + session_name + '/func/sub-' + participant_id + '_' + session_name + '_task-rest_desc-preproc_bold.nii.gz'
    else:
        niftiPath = False

    return niftiPath

def csv2dict(participants_path):
    '''

    returns an array of dictionaries

    #with open(dhcp_dir + 'dhcp_fmri_pipeline/participants.tsv') as csvfile:#, newline='') as csvfile:
    with open(participants_path, newline='\n') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')#, quotechar='')
        for row in reader:
            key = row[0]
            if key in participants:
                print("duplicado")
                pass
    '''
    #participants = list()
    participants = list()
    with open(participants_path) as input_file:
        reader = csv.DictReader(input_file,  delimiter='\t')
        for row in reader:
            #participants.append(row)
            participants.append(row) # table is a numpy array of dictionaries.
    #return np.asarray(participants)
    return participants

def csv2pd(participants_path):
    '''

    :param dictionary:
    :return: returns an array where [subject_index][field_index][0] is field name and [subject_index][field_index][1] is value
    '''

    result = list()

    with open(participants_path, newline='\n') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')  # , quotechar='')
        #headers = next(reader)[0]
        headers = next(reader)
        if len(headers) == 1: # Needed for malformed csv files whose headers can't be read.
            headers = headers[0].split(",")
            headers = list(headers)
            for row in reader:
                line_list = row[0].split(",")
                result.append(line_list)
        else:
            for row in reader:
                result.append(row)
    dataframe = pd.DataFrame(result, columns=headers)
    return dataframe

def get_longeststreak(participant_fd, threshold):
    '''

    :param participant_fd: framewise displacement has to be flaot
    :param threshold: float
    :return: longest streak
    '''
    thresholded_array = (participant_fd['framewise_displacement'] < threshold).values
    i = amount_true = streak = 0
    while i < (len(thresholded_array) -1) :
        if thresholded_array[i]:
            changes = np.diff(thresholded_array[i:])
            amount_true = np.argmax(changes)
            if sum(np.logical_not(changes)) == len(changes):  # Important to try and eliminate the need for this if
                amount_true = len(changes)
            streak = max([streak, amount_true])
            i = i + amount_true
        else:
            i = i + 1
        if amount_true == 0:
            i = i + 1

    return streak

def boldstats(subjects_directory = "/N/project/baby_ICA/dhcp_fmri_cropped", save_pickle = False, save_csv = False):
	'''
	subjects_directory = "/N/project/baby_ICA/dhcp_fmri_cropped" - directory of fmri data 
	save_pickle = False - If not false saves data to pickle in save_pickle
	save_csv = False - If not false saves data to csv in save_csv
	Calls for FSL on linux and returns a pandas dataframe with the stats of the bold timeseries
	 
	
	
	
	'''
	
	#dhcp_root = "/N/project/baby_ICA"
	subject_list = os.listdir(subjects_directory)

	pl = pd.DataFrame({"subject": subject_list})
	pl.insert(1, "mean", -1)
	pl.insert(2, "sd", -1)
	pl.insert(3, "min", -1)
	pl.insert(4, "max", -1)

	for subject in subject_list:
		process = subprocess.Popen(
			["fslstats " + get_niftiPath(subject) + " -m -s -R"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout, stderr = process.communicate()
		outline = stdout.decode("utf-8").split(' ')
		pl.loc[pl["subject"] == subject, "mean"] = float(outline[0])
		pl.loc[pl["subject"] == subject, "sd"] = float(outline[1])
		pl.loc[pl["subject"] == subject, "min"] = float(outline[2])
		pl.loc[pl["subject"] == subject, "max"] = float(outline[3])
	
	if save_pickle:
		pl.to_pickle(save_pickle)
	
	if save_csv:
		pl.to_csv(save_csv)

	return(pl)

def get_scanage(subid):
    '''
        Returns the age at scan for the first session of the subject.
        TODO: handle multiple sessions
        input:
            subid: str - subject id
        output:
            age: float - returns age at scan, 0 if not found.
    '''
    age = 0
    subid = 'sub-' + subid # Incomplete name in subject_table.
    subject_path = dhcp_root() + "/dhcp_fmri_pipeline/" + subid
    if os.path.exists(subject_path):
        for directory in os.listdir(subject_path):
            if directory[0:3] == 'ses':
                amount_sessions = 1 # Change for multiple sessions.
        if amount_sessions < 1:
            print("Session not found.")
        else:
            session_info = subject_path + '/' + subid + '_sessions.tsv'
            with open(session_info, newline='\n') as csvfile:
                reader = csv.reader(csvfile, delimiter='\t')  # , quotechar='')
                headers = next(reader) # Read headers for future functionality
                first_session = next(reader)
                age = first_session[1]
    else:
        print("Subject not found.")

    return float(age)

def get_sesid(subid):
    '''
        Returns the age at scan for the first session of the subject.
        TODO: handle multiple sessions
        input:
            subid: str - subject id
        output:
            age: float - returns age at scan, 0 if not found.
    '''
    age = 0
    subid = 'sub-' + subid # Incomplete name in subject_table.
    subject_path = dhcp_root() + "/dhcp_fmri_pipeline/" + subid
    if os.path.exists(subject_path):
        for directory in os.listdir(subject_path):
            if directory[0:3] == 'ses':
                sesid = directory
                amount_sessions = 1 # Change for multiple sessions.
        if amount_sessions < 1:
            print("Session not found.")
    else:
        print("Subject not found.")

    return sesid

def generate_fullpaths(sufix, subject_list):
    '''
        Generates a txt file with the full paths to the dHCP data type specified
        of the subjects in the pickle file 'subject_list'
        input:
            sufix: str - what comes after subid-sesid.
            subject_list: pkl - pkl file with a list of subjects.
    '''

