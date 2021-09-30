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


#map(folders, make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021.csv") #needs to be run once

#STEP 2: Read in all files in a directory and rename the channels.
#Results in a gating set containing all .fcs files, associated experiment details, and marker details
#Author: Julie

#Do this this for files in timepoint01 directory
#Grace and I decided to write the experiment-markers.csv to the parent directory ie) FSC_files not the timepoint subdirectory.

exp_details_path = list.files(path = paste0(folders[1]), pattern = "_experiment_details.csv", full.names = T) #a way to stay in the parent directory but access the timepoint subdirectories as needed when making gating sets in cyto_setup()

timept01_gating_set <- cyto_setup(path=folders[1], restrict=TRUE, select="fcs", details=F) #details=F; use interactive GUI to paste in experiment details for first timepoint

cyto_details(timept01_gating_set) #currently ONLY the name column

#annotate the experiment details file associated with the gating set using pData.
#our workaround to cyto_setup() not being able to read in the experimental-details.csv file we generated in STEP 1.
tp01_experiment_details <- read_csv(exp_details_path) #import experiment-details.csv
pData(timept01_gating_set)$name<-tp01_experiment_details$name
flowWorkspace::pData(timept01_gating_set)$sample<-tp01_experiment_details$sample
cyto_details(timept01_gating_set) #check what the experiment details for this gs look like
flowWorkspace::pData(timept01_gating_set)$`Outflow Well`<-tp01_experiment_details$`Outflow well`
flowWorkspace::pData(timept01_gating_set)$Media<-tp01_experiment_details$Media
flowWorkspace::pData(timept01_gating_set)$Strain<-tp01_experiment_details$Strain
flowWorkspace::pData(timept01_gating_set)$Type<-tp01_experiment_details$Type
flowWorkspace::pData(timept01_gating_set)$Description<-tp01_experiment_details$Description
flowWorkspace::pData(timept01_gating_set)$generation<-tp01_experiment_details$generation

cyto_details(timept01_gating_set) #all experimental details metadata now associated with the gating set

file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv") #rename the experiment-markers.csv file from cytoexplorer's default to whatever you want. Since this is a universal file to be used across all timepoints I gave it the experiment name. Grace and I decided to have it sit in the parent directory since it's universal file.

#STEP 3:  Perform gating on gating set
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.
#Author: Titir

#First transform the data
timept01_transformed <- cyto_transformer_logicle(timept01_gating_set,
                                              channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A")) # error about channel not being valid can occur if markers = F in cyto_setup() because markers are not assigned
# check if cyto_markers() is NULL
transformed_timept01 <- cyto_transform(timept01_gating_set,
                                       trans = timept01_transformed)


##Gating using the entire timepoint1 dataset
#First we gate for the cells
cyto_gate_draw(transformed_timept01,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = "cytek_gating_all.csv",
)

#Then we define the singlets based on forward scatter height and width
cyto_gate_draw(transformed_timept01,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "cytek_gating_all.csv"
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
               gatingTemplate = "cytek_gating_all.csv",
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
#2 Read in experiment details files using pData
#3 Specify experiment markers
#4 Transform gating set
#5 Apply existing gating file using cyto_gatingTemplate_apply
#6.Output stats file as .csv
#Author: David & Julie

#to be executed from the parent directory

#timept02_gating_set <- cyto_setup(path=folders[2], restrict=TRUE, select="fcs", details=F) #details=F; use interactive GUI to paste in experiment details for first timepoint



analyze_all_exp(folders[-1], markers.new, gating_template
#one folder name,
# we want this function to applied to many folders
#markers and gating template is universal
map(folders[-1], analyze_all_exp(experiment_markers = markers.new, gating_template = "cytek_gating_all.csv"))

#Dev
folder_name <- '02_EE_GAP1_ArchMuts_2021_062121_g21_TD'
#experiment_details <- '02_EE_GAP1_ArchMuts_2021_experiment_details.csv'

my_markers<-c("GFP") #list marker name
channel<-c("B2-A") #list channels
names(my_markers)<-channel

gating_template <- 'cytek_gating_all.csv'

analyze_all_exp = function(folder_name, experiment_markers, gating_template) {
  #we don't need to give experiment_details as one of the arguments in the function, just the folder name is needed and markers, and gating template. from the folder name inside the function it can go one folder down to find the experiment_details.csv files by pattern searching I hope.

  my_path <- paste0("./", folder_name) #gets relative path name for folder to be analyzed

  prefix <- folder_name %>% str_extract("([0-9])+_EE_GAP1_ArchMuts_2021") #extracts the time point number from folder name

  my_expt_details_path <- paste0(my_path,"/",prefix,"_experiment_details.csv") #gets experiment details .csv from correct directory
#  my_experiment_markers <- 'EE_GAP1_ArchMuts_2021-Experiment-Markers.csv' # provided as argument and in working directory

#  my_gating_template <- 'cytek_gating_all.csv' #provided as argument and in working directory

  #  gt <- gatingTemplate("Cytek_gating_all.csv") #redundant code?

  #1. read in files and make a gating set
  timepoint_gating_set <- cyto_setup(path=my_path, select="fcs", details=F, markers = F) #do not set gatingTemplate = gating_template it will overwrite any existing gating template with a BLANK csv file #read in details and markers later

  #2. read in experiment details
  my_experiment_details <- read_csv(my_expt_details_path) #import experiment-details.csv
  flowWorkspace::pData(timepoint_gating_set)$name<-my_experiment_details$name
  flowWorkspace::pData(timepoint_gating_set)$sample<-my_experiment_details$sample
  flowWorkspace::pData(timepoint_gating_set)$`Outflow Well`<-my_experiment_details$`Outflow well`
  flowWorkspace::pData(timepoint_gating_set)$Media<-my_experiment_details$Media
  flowWorkspace::pData(timepoint_gating_set)$Strain<-my_experiment_details$Strain
  flowWorkspace::pData(timepoint_gating_set)$Type<-my_experiment_details$Type
  flowWorkspace::pData(timepoint_gating_set)$Description<-my_experiment_details$Description
  flowWorkspace::pData(timepoint_gating_set)$generation<-my_experiment_details$generation

  #3. specify markers
  #cyto_markers(timepoint_gating_set) #currently no markers
  markernames(timepoint_gating_set)<-my_markers
  #cyto_markers(timepoint_gating_set) GFP channel now specified

  #4. transform data
  timepoint_gating_set_transformed <- cyto_transformer_logicle(timepoint_gating_set,
                                                   channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A")) #transforms but returns the gating set as a list
  transformed_timepoint_gating_set <- cyto_transform(timepoint_gating_set,
                                         trans = timepoint_gating_set_transformed) # applies the transformation and converts the list to a GatingSet object

  #apply gating-template.csv to transformed gating set
  cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate = gating_template)

  cyto_plot_gating_scheme(transformed_timepoint_gating_set)[[,"Strain"=="DGY500"]# plot with drawn gates will appear if gates were indeed applied
  cyto_plot_gating_scheme(transformed_timepoint_gating_set)["Strain"=="DGY1"]
  cyto_plot_gating_scheme(transformed_timepoint_gating_set)["Strain"=="DGY1315"]
  #, #DGY500 sample
   #                       gatingTemplate = "time01_gating.csv",
    #                      back_gate = TRUE,
     #                     gate_track = TRUE)

  #write stats freq file for % of cells inside each gate
  stats_timept2 <- cyto_stats_compute(transformed_timepoint_gating_set,
                                      parent = "Single_cells",
                                      alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                      stat="freq",
                                      #gate = gt,
                                      save_as = paste0("stats_timept_",prefix,".csv") #writes to working directory
                                      )
  }

#example

analyze_all_exp('2_EE_GAP1_ArchMuts_2021_062121_g21_TD',
        #        '2_EE_GAP1_ArchMuts_2021_experiment_details.csv',
         #       'EE_GAP1_ArchMuts_2021-Experiment-Markers.csv',
                'cytek_gating_all.csv')

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Grace

map(folders[-1], analyze_all_exp(experiment_markers = markers.new, gating_template = "cytek_gating_all.csv")) #sorry Grace I got ahead of myself. Does

#STEP 7:  Combine stats_freq.csv files into a single dataframe
#Pull in all stats_freq files from directories and assemble into a single dataframe
#Author: Julie

#STEP 8: Assess gates
#Determine whether =>95% of controls are in the correct gate
#Author: Julie





