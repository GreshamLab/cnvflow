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

# Set working directory and get list of subdirectories
#setwd('./Summer 2021 Group LTEE/FCS files/') #Grace's working directory
#setwd('/Volumes/GoogleDrive/My Drive/Gresham Lab_Papers/2021/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files') #David's working directory
#setwd('G:/.shortcut-targets-by-id/1Bioj1YP_I7P8tqgmg4Zbt4EAfhb7J-0w/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files') #Titir's working directory

folders = list.dirs()[-1]

#STEP 1: Generate experiment details file.
#A .csv file that contains the list of .fcs files in the directory and the associated metadata for each sample
#Author: Grace

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

map(folders, make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021.csv") #needs to be run once

#STEP 2: Read in all files in a directory and rename the channels.
#Results in a gating set containing all .fcs files, associated experiment details, and marker details
#Author: Julie

exp_details_path = list.files(path = paste0(folders[1]), pattern = "_experiment_details.csv", full.names = T)

timept01_gating_set <- cyto_setup(path=folders[1], restrict=TRUE, select="fcs", details=F)

#use pData to annotate the experiment details file associated with the gating set
tp01_experiment_details <- read_csv(exp_details_path) #import experiment-details.csv
pData(timept01_gating_set)$name<-tp01_experiment_details$name
flowWorkspace::pData(timept01_gating_set)$sample<-tp01_experiment_details$sample
flowWorkspace::pData(timept01_gating_set)$`Outflow Well`<-tp01_experiment_details$`Outflow well`
flowWorkspace::pData(timept01_gating_set)$Media<-tp01_experiment_details$Media
flowWorkspace::pData(timept01_gating_set)$Strain<-tp01_experiment_details$Strain
flowWorkspace::pData(timept01_gating_set)$Type<-tp01_experiment_details$Type
flowWorkspace::pData(timept01_gating_set)$Description<-tp01_experiment_details$Description
flowWorkspace::pData(timept01_gating_set)$generation<-tp01_experiment_details$generation

file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv") #rename the experiment-markers.csv file. Need to do once.

#STEP 3:  Perform gating on gating set
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.
#Author: Titir

#Log transform the data
timept01_transformed <- cyto_transformer_log(timept01_gating_set,
                                                 channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A")) #returns it as a list
transformed_timept01 <- cyto_transform(timept01_gating_set,
                                       trans = timept01_transformed) #applies the the transformation and returns it as a gatingSet

##Gating using the entire timepoint01 dataset.
#First we gate for the cells
cyto_gate_draw(transformed_timept01,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = "cytek_gating.csv",
)

#Then we define the singlets based on forward scatter height and width
cyto_gate_draw(transformed_timept01,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "cytek_gating.csv"
)

#Gating for CNVs using the 0,1 and 2 copy controls:
zero_copy <- cyto_extract(transformed_timept01, "Single_cells")[[30]] #DGY1

one_copy <- cyto_extract(transformed_timept01, "Single_cells")[[1]] #DGY500

two_copy <- cyto_extract(transformed_timept01, "Single_cells")[[31]] #DGY1315

cyto_gate_draw(transformed_timept01,
               parent = "Single_cells",
               alias = c("zero_copy", "one_copy", "two_copy","multi_copy"), #defines gate names
               channels = c("FSC-A","B2-A"),
               axes_limits = "data",
#               select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
               gatingTemplate = "cytek_gating.csv",
               overlay = c(zero_copy, one_copy, two_copy),
               point_col = c("black", "green", "red", "blue")
)

#STEP 4:  Generate statistics tables
#Results in a .csv file in tidy format that includes all metadata and specifies proportion of cells with 0, 1, 2, 3+ copies, a .csv file with overall median GFP and FSC-A, and a third .csv file with gate-wise median GFP and FSC-A.
#Author: Titir

stats_freq_01 <- cyto_stats_compute(transformed_timept01,
                                  parent = "Single_cells",
                                  alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                  stat="freq",
                                  save_as = "stats_freq_01.csv")

stats_median_overall_01 <- cyto_stats_compute(transformed_timept01,
                                      parent = c("Single_cells"),
                                      alias = c("Single_cells"),
                                      channels = c("FSC-A", "B2-A"),
                                      stat="median",
                                      save_as = "stats_median_01_overall.csv")

stats_median_gatewise_01 <- cyto_stats_compute(transformed_timept01,
                                              parent = c("Single_cells"),
                                              alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                              channels = c("FSC-A", "B2-A"),
                                              stat="median",
                                              save_as = "stats_median_01_gatewise.csv")

#STEP 5:  Use function to perform analysis
#A function that will
#1 Read in all the files in a folder
#2 Read in experiment details files using pData
#3 Specify experiment markers
#4 Transform gating set
#5 Apply existing gating file using cyto_gatingTemplate_apply
#6.Output stats file as .csv
#Author: David & Julie

my_markers<-c("GFP") #list your marker name(s)
channel<-c("B2-A") #list your channel(s)
names(my_markers)<-channel

analyze_all_exp = function(folder_name, my_markers,"cytek_gating.csv") {

  my_path <- folder_name #gets relative path name for folder to be analyzed

  prefix <- folder_name %>% str_extract("([0-9])+_EE_GAP1_ArchMuts_2021") #extracts the time point number from folder name

  my_expt_details_path <- paste0(my_path,"/",prefix,"_experiment_details.csv") #gets experiment details .csv from correct directory

  #1. read in files and make a gating set
  timepoint_gating_set <- cyto_setup(path=my_path, select="fcs", details=F, markers = F)

  #2. read in experiment details for that gating set
  my_experiment_details <- read_csv(my_expt_details_path) #import experiment-details.csv
  flowWorkspace::pData(timepoint_gating_set)$name<-my_experiment_details$name
  flowWorkspace::pData(timepoint_gating_set)$sample<-my_experiment_details$sample
  flowWorkspace::pData(timepoint_gating_set)$`Outflow Well`<-my_experiment_details$`Outflow well`
  flowWorkspace::pData(timepoint_gating_set)$Media<-my_experiment_details$Media
  flowWorkspace::pData(timepoint_gating_set)$Strain<-my_experiment_details$Strain
  flowWorkspace::pData(timepoint_gating_set)$Type<-my_experiment_details$Type
  flowWorkspace::pData(timepoint_gating_set)$Description<-my_experiment_details$Description
  flowWorkspace::pData(timepoint_gating_set)$generation<-my_experiment_details$generation

  #3. specify markers for that gating set
  markernames(timepoint_gating_set)<-my_markers

  #4. transform data
  timepoint_gating_set_transformed <- cyto_transformer_log(timepoint_gating_set,
                                                           channels =c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A")) #transforms but returns the gating set as a list
  transformed_timepoint_gating_set<- cyto_transform(timepoint_gating_set,
                                                    trans = timepoint_gating_set_transformed)

  #5. apply gating-template.csv to transformed gating set
  cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= 'cytek_gating.csv')

  #6. write stats: freq file for % of cells inside each gate, median FSC and GFP for each population, median FSC and GFP for each gated population
  #Titir
  stats_freq <- cyto_stats_compute(transformed_timepoint_gating_set,
                                      parent = c("Single_cells"),
                                      alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                      stat="freq",
                                      save_as = paste0("stats_freq_",prefix,".csv") #writes to working directory
                                      )
  stats_median_overall <- cyto_stats_compute(transformed_timepoint_gating_set,
                                     parent = c("Single_cells"),
                                     alias  = c("Single_cells"),
                                     channels = c("FSC-A", "B2-A"),
                                     stat="median",
                                     save_as = paste0("stats_median_overall", prefix,".csv"))

  stats_median_gatewise <- cyto_stats_compute(transformed_timepoint_gating_set,
                                              parent = c("Single_cells"),
                                              alias  = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                              channels = c("FSC-A", "B2-A"),
                                              stat="median",
                                              save_as = paste0("stats_median_gatewise", prefix,".csv"))
}

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Julie

map(folders[-1], analyze_all_exp, my_markers, "cytek_gating.csv")

#STEP 7:  Combine stats_freq.csv and stats_median.csv files into a single dataframe
#Pull in all stats_* files from directories and assemble into a single dataframe
#Author: Julie

list.files(path = ".", pattern = "stats_freq") %>%
  read_csv() %>%
  write_csv(file = "stats_freq_all_timepoints.csv")

list.files(path = ".", pattern = "stats_median_overall") %>%
  read_csv() %>%
  write_csv(file = "stats_median_overall_all_timepoints.csv")

list.files(path = ".", pattern = "stats_median_gatewise") %>%
  read.csv() %>%
  write_csv(file = "stats_median_gatewise_all_timepoints.csv")

#STEP 8: Plot time series & assess gates
#Determine whether =>90-95% of controls are in the correct gate
#Author: Grace

#dplyr::select filter controls, zero copy, one copy, two copy,






