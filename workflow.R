#workflow.R
#Started: Sept 21, 2021
#Authors: Julie Chuong, Titir De, Grace Avecilla, David Gresham

#STEP 1: Generate experiment details file.
#A .csv file that contains the list of .fcs files in the directory and the associated metadata for each sample
#Author: Grace


#STEP 2: Read in all files in a directory and rename the channels.
#Results in a gating set containing all .fcs files, associated experiment details, and marker details
#Author: Julie


#STEP 3:  Perform gating on gating set
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.
#Author: Titir


#STEP 4:  Generate statistics table
#Results in a .csv file in tidy format that includes all metadata and specifies proportion of cells with 0, 1, 2, 3+ copies]
#Author: Titir

#STEP 5:  Use function to perform analysis
#A function that will
#1 Read in all the files in a folder
#2 Red in experiment details files
#3 Specify experiment markers
#4 Transform gating set
#5 Apply existing gating file
#6.Output stats file as .csv
#Author: David

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Grace

#STEP 7:  Combine stats_freq.csv files into a single dataframe
#Pull in all stats_freq files from directories and assemble into a single dataframe
#Author: Julie

#STEP 8: Assess gates
#Determine whether =>95% of controls are in the correct gate
#Author: Julie





