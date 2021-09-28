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

setwd('/Volumes/GoogleDrive/My Drive/Gresham Lab_Papers/2021/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files')

exp_details_path = list.files(path = paste0(folders[1]), pattern = "_experiment_details.csv", full.names = T) #a way to stay in the parent directory but access the timepoint subdirectories as needed when making gating sets in cyto_setup()

timept01_gating_set <- cyto_setup(path=folders[1], restrict=TRUE, select="fcs", details=F) #details=F; use interactive GUI to paste in experiment details for first timepoint

#annotarte using pData.  Figure out how to add the whole dataframe
test_1_EE_GAP1_ArchMuts_2021_experiment_details <- read_csv("01_EE_GAP1_ArchMuts_2021_061621_g8_GA/1_EE_GAP1_ArchMuts_2021_experiment_details.csv")
pData(transformed_timept01)$sample <- test_1_EE_GAP1_ArchMuts_2021_experiment_details$sample

file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv") #rename the experiment-markers.csv file from cytoexplorer's default to whatever you want. Since this is a universal file to be used across all timepoints I gave it the experiment name. Grace and I decided to have it sit in the parent directory since it's universal file.

#STEP 3:  Perform gating on gating set
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.
#Author: Titir

#First transform the data
timept01_transformed <- cyto_transformer_logicle(timept01_gating_set,
                                              channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A"))
transformed_timept01 <- cyto_transform(timept01_gating_set,
                                       trans = timept01_transformed)

##Gating using the entire timepoint1 dataset
#First we gate for the cells
cyto_gate_draw(transformed_timept01,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = "Cytek_gating.csv",
)

#Then we define the singlets based on forward scatter height and width
cyto_gate_draw(transformed_timept01,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "Cytek_gating.csv"
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
               gatingTemplate = "Cytek_gating.csv",
               overlay = c(zero_copy, one_copy, two_copy),
               point_col = c("black", "green", "red", "blue")
)

#STEP 4:  Generate statistics table
#Results in a .csv file in tidy format that includes all metadata and specifies proportion of cells with 0, 1, 2, 3+ copies]
#Author: Titir

stats_timept1 <- cyto_stats_compute(transformed_timept01,
                                  parent = "Single_cells",
                                  alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                  stat="freq",
                                  save_as = "stats_timept1.csv")

#STEP 5:  Use function to perform analysis
#A function that will
#1 Read in all the files in a folder
#2 Read in experiment details files
#3 Specify experiment markers
#4 Transform gating set
#5 Apply existing gating file
#6.Output stats file as .csv
#Author: David

#to be executed from the parent directory

timept02_gating_set <- cyto_setup(path=folders[2], restrict=TRUE, markers=F, select="fcs", details=F) #details=F; use interactive GUI to paste in experiment details for first timepoint


analyze_all_exp = function(folder_name, experiment_details, experiment_markers, gating_template) {

  #Dev
  folder_name <- '2_EE_GAP1_ArchMuts_2021_062121_g21_TD'
  experiment_details <- '2_EE_GAP1_ArchMuts_2021_experiment_details.csv'
  experiment_markers <- 'EE_GAP1_ArchMuts_2021-Experiment-Markers.csv'
  gating_template <- 'Cytek_gating.csv'

  my_path <- paste0("./", folder_name) #gets relative path name for folder to be analyzed

  prefix <- folder_name %>% str_extract("([0-9])+_EE_GAP1_ArchMuts_2021") #extracts the time point number from folder name

  my_expt_details <- paste0(my_path,"/",experiment_details) #gets experiment details .csv from correct directory

  gt <- gatingTemplate("Cytek_gating.csv")
#  my_gates <- cyto_gate_extract("Cells", gatingTemplate="Cytek_gating.csv")

#  my_experiment_markers <- 'EE_GAP1_ArchMuts_2021-Experiment-Markers.csv' # provided as argument and in working directory

#  my_gating_template <- 'Cytek_gating.csv' #provided as argument and in working directory

  timepoint_gating_set <- cyto_setup(path=my_path, restrict=TRUE, select="fcs", details=F, gatingTemplate = gating_template)

  #transform data

  timepoint_gating_set_transformed <- cyto_transformer_logicle(timepoint_gating_set,
                                                   channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A"))
  transformed_timepoint_gating_set <- cyto_transform(timepoint_gating_set,
                                         trans = timepoint_gating_set_transformed) #

  stats_timept2 <- cyto_stats_compute(transformed_timepoint_gating_set,
                                      parent = "Single_cells",
                                      alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                      stat="freq",
                                      gate = gt,
                                      save_as = paste0("stats_timept",prefix,".csv") #writes to working directory
                                      )
  }

#example

analyze_all_exp('2_EE_GAP1_ArchMuts_2021_062121_g21_TD',
                '2_EE_GAP1_ArchMuts_2021_experiment_details.csv',
                'EE_GAP1_ArchMuts_2021-Experiment-Markers.csv',
                'Cytek_gating.csv')

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Grace

#STEP 7:  Combine stats_freq.csv files into a single dataframe
#Pull in all stats_freq files from directories and assemble into a single dataframe
#Author: Julie

#STEP 8: Assess gates
#Determine whether =>95% of controls are in the correct gate
#Author: Julie





