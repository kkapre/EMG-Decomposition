# EMG-Decomposition
Decomposing raw electromyography data into motor unit action potentials (MUAPs) can help with understanding the science of muscles and with clinical diagnoses.  

Based off project for Neuroengineerig class at UCLA

### EMG_Background.pdf
Gives background on EMG data and EMG decomposition

### EMG_Decomposition.m
MATLAB code starts with raw EMG data and outputs the shape of spikes from individual motor units through the process of filtering, detecting and aligning spikes, and k-means clustering. It also gives firing rates of each detected motor unit. 

Place data "EMG_example_2_fs_2k.csv" in the same directory to run. 
### EMG_Decomposition_Report.pdf
Analysis sumamry with figures and results
