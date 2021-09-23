#workflow.R
#Started: Sept 21, 2021
#Authors: Julie Chuong, Titir De, Grace Avecilla, David Gresham


##Install CytoExplorer package and requirements (can be skipped if already installed)
#library(BiocManager)
#install("cytolib", "flowCore", "flowWorkspace", "openCyto")

##Install CytoExploreR from GitHub:
#library(devtools)
#devtools::install_github("DillonHammill/CytoExploreR")

# Load required packages
library(CytoExploreR)
library(tidyverse)

#STEP 1: Generate experiment details file.
#A .csv file that contains the list of .fcs files in the directory and the associated metadata for each sample
#Author: Grace
setwd('./Summer 2021 Group LTEE/FCS files/')

setwd('./Summer 2021 Group LTEE/FCS files/')

make_exp_details = function(folder_name, samplesheet) {
  pref = folder_name %>% str_extract("([0-9])+_EE_GAP1_ArchMuts_2021")
  generation = folder_name %>% str_extract("[g]\\d+") %>% str_remove("g")

  files = as_tibble(list.files(paste0(folder_name))) %>%
    separate(value, into = c("well", "samp"), sep = " ", remove = F) %>%
    mutate(well = str_extract(well, "([A-Z])([0-9]){1,2}$")) %>%
    mutate(samp = str_remove(samp, ".fcs")) %>%
    mutate(sample = case_when(str_detect(value, "Unstained") ~ "ctrl0",
                              str_detect(value, "DGY500") ~ "ctrl1",
                              str_detect(value, "DGY1315") ~ "ctrl2",
                              TRUE ~ samp)) %>%
    select(value,sample) %>%
    rename(name = value) %>%
    filter(!is.na(sample))

  all = files %>%
    left_join(read_csv(paste0("./",samplesheet)), by = c("sample" = "Sample name")) %>%
    mutate(generation = as.numeric(generation))

  write_csv(all, file = paste0(folder_name,"/",pref,"_experiment_details.csv"))

}

folders = list.dirs()[-1]
map(folders, make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021.csv")

#STEP 2: Read in all files in a directory and rename the channels.
#Results in a gating set containing all .fcs files, associated experiment details, and marker details
#Author: Julie

#Do this this for files in timepoint01 directory
#Grace and I decided to write the experiment-markers.csv to the parent directory ie) FSC_files not the timepoint subdirectory.

exp_details_path = list.files(path = paste0(folders[1]), pattern = "_experiment_details.csv", full.names = T) #a way to stay in the parent directory but access the timepoint subdirectories as needed when making gating sets in cyto_setup()

timept01_gating_set <- cyto_setup(path=folders[1],restrict=TRUE, select="fcs", details=F) #details=F because experiment-details.csv files were already generated in STEP 1

file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv") #rename the experiment-markers.csv file from cytoexplorer's default to whatever you want. Since this is a universal file to be used across all timepoints I gave it the experiment name. Grace and I decided to have it sit in the parent directory since it's universal file.

#STEP 3:  Perform gating on gating set
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.
#Author: Titir

#First we transform the data
timept01_transformed <- cyto_transformer_logicle(timept01_gating_set,
                                              channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A"))


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





