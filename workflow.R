#workflow.R
#Started: Sept 21, 2021
#Authors: Julie Chuong, Titir De, Grace Avecilla, David Gresham

##Install CytoExplorer package and requirements (can be skipped if already installed)
#library(BiocManager)
#install.packages("cytolib", "flowCore", "flowWorkspace", "openCyto")

##Install CytoExploreR from GitHub:
#library(devtools)
#devtools::install_github("DillonHammill/CytoExploreR")

# Load required packages
library(CytoExploreR)
library(tidyverse)
library(purrr)
library(ggridges)

# Set working directory and get list of subdirectories
#setwd('../FCS files/') #Grace's working directory
#setwd('/Volumes/GoogleDrive/My Drive/Gresham Lab_Papers/2021/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files') #David's working directory
#setwd('G:/.shortcut-targets-by-id/1Bioj1YP_I7P8tqgmg4Zbt4EAfhb7J-0w/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files') #Titir's working directory
#setwd("/Volumes/GoogleDrive/My Drive/greshamlab/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files") #Julie's WD
setwd("/Volumes/GoogleDrive/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")  #Julie's WD
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
#map(folders[1], make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021.csv")
#STEP 2: Read in all files in a directory and rename the channels.
#Results in one timepoint's gating set containing all .fcs files, associated experiment details, and marker details
#Author: Julie
folders = list.dirs()[-1]
exp_details_path = list.files(path = paste0(folders[4]), pattern = "_experiment_details.csv", full.names = T)
#exp_details_path = list.files(path = folder08, pattern = "_experiment_details.csv", full.names = T)
timepoint_gating_set <- cyto_setup(path = paste0(folders[4]), restrict=TRUE, select="fcs", details=F) #edit Markers on Viewer pane, Save & Close

#use pData to annotate the experiment details file associated with the gating set
experiment_details <- read_csv(exp_details_path) #import experiment-details.csv
for(i in 1:length(names(experiment_details))){
  flowWorkspace::pData(timepoint_gating_set)[names(experiment_details[i])]<-experiment_details[i]
  }
#flowWorkspace::pData(timepoint_gating_set)$name<-experiment_details$name
#flowWorkspace::pData(timepoint_gating_set)$sample<-experiment_details$sample
#flowWorkspace::pData(timepoint_gating_set)$`Outflow Well`<-experiment_details$`Outflow well`
#flowWorkspace::pData(timepoint_gating_set)$Media<-experiment_details$Media
#flowWorkspace::pData(timepoint_gating_set)$Strain<-experiment_details$Strain
#flowWorkspace::pData(timepoint_gating_set)$Type<-experiment_details$Type
#flowWorkspace::pData(timepoint_gating_set)$Description<-experiment_details$Description
#flowWorkspace::pData(timepoint_gating_set)$generation<-experiment_details$generation

#file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv") #rename the experiment-markers.csv file. Need to do once.

#STEP 3:  Perform gating on gating set
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.
#Author: Titir

#Log transform the data
# looks useful if I want to choose different transformation: https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR-Transformations.html
#timept_transformed <- cyto_transformer_log(timepoint_gating_set,
#                      channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A")) #returns it as a list
#transformed_timepoint_gating_set <- cyto_transform(timepoint_gating_set,
#                      trans = timept_transformed) #applies the the transformation and returns it as a gatingSet
#timept_transformed <- cyto_transformer_logicle(timepoint_gating_set,
#                                           channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A"),
#                                           widthBasis = -10) #returns it as a list
#transformed_timepoint_gating_set <- cyto_transform(timepoint_gating_set,
#                                     trans = timept_transformed) #applies the the transformation and returns it as a gatingSet
GFP_trans <- cyto_transformer_logicle(timepoint_gating_set,
                                      channels = c("B2-A"),
                                      widthBasis = -10
)#returns it as a list
FSC_SSC_trans <- cyto_transformer_log(timepoint_gating_set,
                                      channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H")
)
combined_trans <- cyto_transformer_combine(GFP_trans,FSC_SSC_trans)
transformed_timepoint_gating_set <- cyto_transform(timepoint_gating_set,
                                                   trans = combined_trans) #applies the the transformation and returns it as a gatingSet


#quickly check the transformation by plotting the data
#cyto_plot_explore(transformed_timepoint_gating_set[c(2,14,16,17,18,19,21)],
#                  channels_x = "FSC-A",
#                  channels_y = "GFP",
#                  axes_limits = "data")

##Gating using the entire timepoint dataset.

#if you already have a gating template and don't need to draw gates, then skip drawing and apply the gating template to your gating set
cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= "cytek_gating_01_02_04_v2.csv")

#First we gate for the cells
cyto_gate_draw(transformed_timepoint_gating_set,
  #transformed_logicle_timept,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = "cytek_gating_01_02_04_v2.csv",
)

#Then we define the singlets based on forward scatter height and width
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "cytek_gating_01_02_04_v2.csv"
)

#Gating for CNVs using the 0,1 and 2 copy controls:
DGY1 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(30,61,92)] #DGY1

DGY500 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(1,32,63)] #DGY500

DGY1315 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(31,62,93)] #DGY1315

cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Single_cells", #first color
               alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
               channels = c("FSC-A","B2-A"),
               axes_limits = "data",
               #select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
               gatingTemplate = "cytek_gating_01_02_04_v2.csv",
               overlay = c(DGY1, DGY500, DGY1315),
               point_col = c("gray", "green", "red", "blue")
)

#STEP 4:  Generate statistics tables
#Results in a .csv file in tidy format that includes all metadata and specifies proportion of cells with 0, 1, 2, 3+ copies, a .csv file with overall median GFP and FSC-A, and a third .csv file with gate-wise median GFP and FSC-A.
#Author: Titir & Julie

stats_freq_01 <- cyto_stats_compute(transformed_timepoint_gating_set01,
                                  parent = "Single_cells",
                                  alias = c("zero_copy", "one_copy", "two_or_more_copy"),
                                  stat="freq",
                                  save_as = "stats_freq_01.csv")

stats_median_overall_01 <- cyto_stats_compute(transformed_timepoint_gating_set01,
                                      parent = c("Single_cells"),
                                      alias = c("Single_cells"),
                                      channels = c("FSC-A", "B2-A"),
                                      stat="median",
                                      save_as = "stats_median_overall_01.csv")

stats_median_gatewise_01 <- cyto_stats_compute(transformed_timepoint_gating_set01,
                                              parent = c("Single_cells"),
                                              alias = c("zero_copy", "one_copy", "two_or_more_copy"),
                                              channels = c("FSC-A", "B2-A"),
                                              stat="median",
                                              save_as = "stats_median_gatewise_01.csv")

timepoint_raw_list <- cyto_extract(transformed_timepoint_gating_set, parent = "Single_cells", raw = TRUE, channels = c("FSC-A", "B2-A")) #raw flow data of each single cell as a list of matrices

map_df(timepoint_raw_list, ~as.data.frame(.x), .id="name") %>% #convert to df, put list name in new column
  mutate(name = as.factor(name)) %>% #convert `name` to factor
  left_join(experiment_details %>% #join by name column to add metadata
  mutate(generation = as.factor(unique(experiment_details$generation)))) %>%
  mutate(B2A_FSC = `B2-A`/`FSC-A`) %>% #compute normalized fluor
  write_csv(paste0("01_02_04_v2_SingleCellDistributions_",prefix,".csv"))

###################################################
# plot ridgeplots (histograms):
# normalized fluorescence histograms of controls at each generations like in Lauer et al. Fig 2A.
# instead of looking at the MEDIAN GFP values per population per generation,
# we want to see the ENTIRE DISTRIBUTION of GFP of the cells per population per generation.

sc_distributions_g8 <- read.csv("01_02_04_v2_SingleCellDistributions_01_EE_GAP1_ArchMuts_2021.csv", stringsAsFactors = T) %>% mutate(generation = factor(generation, levels = unique(sc_distributions_g8$generation))) %>%
  mutate(name = factor(name, levels = unique(rev(c("Experiment_042-Plate_001-Reference Group-B3 Unstained (Cells).fcs",
                            "Experiment_042-Plate_001-1 copy control-D3 DGY500.fcs",
                            "Experiment_042-Plate_001-Reference Group-F3 DGY1315 mCitrine (Cells).fcs",
                            "Experiment_042-Plate_001-Experimental-H3 gap1_1.fcs",
                            "Experiment_042-Plate_001-Experimental-G4 gap1_2.fcs",
                            "Experiment_042-Plate_001-Experimental-H5 gap1_3.fcs",
                            "Experiment_042-Plate_001-Experimental-G6 gap1_4.fcs",
                            "Experiment_042-Plate_001-Experimental-H7 gap1_5.fcs",
                            "Experiment_042-Plate_001-Experimental-C4 gap1_ltr_1.fcs",
                            "Experiment_042-Plate_001-Experimental-D5 gap1_ltr_2.fcs",
                            "Experiment_042-Plate_001-Experimental-C6 gap1_ltr_3.fcs",
                            "Experiment_042-Plate_001-Experimental-D7 gap1_ltr_4.fcs",
                            "Experiment_042-Plate_001-Experimental-C8 gap1_ltr_5.fcs",
                            "Experiment_042-Plate_001-Experimental-B9 gap1_ltr_6.fcs",
                            "Experiment_042-Plate_001-Experimental-H9 gap1_ltr_7.fcs",
                            "Experiment_042-Plate_001-Experimental-E10 gap1_ltr_8.fcs",
                            "Experiment_042-Plate_001-Experimental-D9 gap1_ars_6.fcs",
                            "Experiment_042-Plate_001-Experimental-E4 gap1_ars_1.fcs",
                            "Experiment_042-Plate_001-Experimental-E6 gap1_ars_3.fcs",
                            "Experiment_042-Plate_001-Experimental-F7 gap1_ars_4.fcs",
                            "Experiment_042-Plate_001-Experimental-E8 gap1_ars_5.fcs",
                            "Experiment_042-Plate_001-Experimental-D9 gap1_ars_6.fcs",
                            "Experiment_042-Plate_001-Experimental-A10 gap1_ars_7.fcs",
                            "Experiment_042-Plate_001-Experimental-G10 gap1_ars_8.fcs",
                            "Experiment_042-Plate_001-Experimental-A4 gap1_all_1.fcs",
                            "Experiment_042-Plate_001-Experimental-A6 gap1_all_3.fcs",
                            "Experiment_042-Plate_001-Experimental-B5 gap1_all_2.fcs",
                            "Experiment_042-Plate_001-Experimental-B7 gap1_all_4.fcs",
                            "Experiment_042-Plate_001-Experimental-A8 gap1_all_5.fcs",
                            "Experiment_042-Plate_001-Experimental-G8 gap1_all_6.fcs",
                            "Experiment_042-Plate_001-Experimental-F9 gap1_all_7.fcs",
                            "Experiment_042-Plate_001-Experimental-C10 gap1_all_8.fcs"
)))))%>%
ggplot(aes(x = B2A_FSC, y = name, fill = Description)) +
  geom_density_ridges(scale=1.5, quantile_lines = TRUE, quantiles = 2) +
  xlab("mCitrine fluorescence/forward scatter") +
  ylab("sample") +
  ggtitle("generation 8 ridgeplots") +
  theme_minimal() +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.5))) + #expands the graph space or else the top is cut off
  scale_fill_discrete(breaks=c("0 copy control",
                               "1 copy control",
                               "2 copy control",
                               "GAP1 WT architecture",
                               "GAP1 LTR KO",
                               "GAP1 ARS KO",
                               "GAP1 LTR + ARS KO"))+ #change order of legend items
  theme(
    legend.text = element_text(family="Arial", size = 12),#edit legend text font and size
    legend.title = element_blank(), #remove legend title
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
    )
ggsave()
#Ridgeplot - combine populations/replicates
ggplot(aes(x = B2A_FSC, y = Description, fill = Description)) +
  geom_density_ridges(scale=1.5) +
  xlab("mCitrine fluorescence/forward scatter") +
  ylab("sample") +
  ggtitle("generation 8 ridgeplots") +
  scale_fill_discrete(breaks=c("0 copy control", #change order of legend items
                               "1 copy control",
                               "2 copy control",
                               "GAP1 WT architecture",
                               "GAP1 LTR KO",
                               "GAP1 ARS KO",
                               "GAP1 LTR + ARS KO")) +
  theme(axis.text.x = element_text(family="Arial", size = 10, color = "black"),
        axis.text.y = element_text(family="Arial", size = 10, color = "black"))
  theme_minimal()




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

analyze_all_exp = function(folder_name, my_markers, gating_template="cytek_gating.csv") {

  path <- folder_name #gets relative path name for folder to be analyzed

  prefix <- folder_name %>% str_extract("([0-9])+_EE_GAP1_ArchMuts_2021") #extracts the time point number from folder name

  exp_details_path <- paste0(path,"/",prefix,"_experiment_details.csv") #gets experiment details .csv from correct directory
  #exp_details_path <- list.files(path = paste0(path), pattern = "_experiment_details.csv", full.names = T)
  #1. read in files and make a gating set
  print(path)
  timepoint_gating_set <- cyto_setup(path=path, select="fcs", details=F, markers = F)

  #2. read in experiment details for that gating set
  experiment_details <- read_csv(exp_details_path, show_col_types = F) #import experiment-details.csv
  #Write For Loop: for column in exp_details_path, add that column to timepoint_gating_set's metadata
  experiment_details <- read_csv(exp_details_path) #import experiment-details.csv
  for(i in 1:length(names(experiment_details))){
    flowWorkspace::pData(timepoint_gating_set)[names(experiment_details[i])]<-experiment_details[i]
  }

  #3. specify markers for that gating set
  markernames(timepoint_gating_set)<-my_markers

  #4. transform data
#  timepoint_gating_set_transformed <- cyto_transformer_log(timepoint_gating_set,
#                                                           channels =c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A")) #transforms but returns the gating set as a list
#  transformed_timepoint_gating_set<- cyto_transform(timepoint_gating_set,
#                                                    trans = timepoint_gating_set_transformed)
  GFP_trans <- cyto_transformer_logicle(timepoint_gating_set,
                                                 channels = c("B2-A"),
                                                 widthBasis = -10
                                                 )#returns it as a list
  FSC_SSC_trans <- cyto_transformer_log(timepoint_gating_set,
                        channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H")
                        )
  combined_trans <- cyto_transformer_combine(GFP_trans,FSC_SSC_trans)
  transformed_timepoint_gating_set <- cyto_transform(timepoint_gating_set,
                                       trans = combined_trans) #applies the the transformation and returns it as a gatingSet
  #cyto_plot_explore(transformed_timepoint_gating_set, #quick plot check that it looks okay
  #                  channels_x = "FSC-A",
  #                  channels_y = "GFP",
  #                  axes_limits = "data")
  #First we gate for the cells
  #cyto_gate_draw(transformed_timepoint_gating_set,
  #               parent = "root",
  #               alias = "Cells",
  #               channels = c("FSC-A","SSC-A"),
  #               axes_limits = "data",
  #               gatingTemplate = "cytek_gating_01_02_04.csv",
  #)

  #Then we define the singlets based on forward scatter height and width
  #cyto_gate_draw(transformed_timepoint_gating_set,
  #               parent = "Cells",
  #               alias = "Single_cells",
  #               channels = c("FSC-A","FSC-H"),
  #               axes_limits = "data",
  #               gatingTemplate = "cytek_gating_01_02_04.csv"
  #)
#  zero_copy <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(30,61)] #DGY1
#  one_copy <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(1,32)] #DGY500
#  two_or_more_copy <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(31,62)] #DGY1315
#  cyto_gate_draw(transformed_timepoint_gating_set,
#                 parent = "Single_cells", #first color
#                 alias = c("zero_copy", "one_copy", "two_or_more_copy","multi_copy"), #defines gate names
#                 channels = c("FSC-A","B2-A"),
#                 axes_limits = "data",
#                 select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
#                 gatingTemplate = "cytek_gating_01_02_04.csv",
#                 overlay = c(zero_copy, one_copy, two_or_more_copy),
#                 point_col = c("black", "green", "red", "blue")
#  )

  #5. apply gating-template.csv to transformed gating set
  cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= gating_template)
#  cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= "cytek_gating_01_02_04_v2.csv")

# cyto_plot_profile(transformed_timepoint_gating_set[1:3],
#                    parent = "Single_cells",
#                    channels = c("FSC-A","GFP"),
#                    legend = TRUE,
#                    legend_text = 4)

  #6. write stats: freq file for % of cells inside each gate, median FSC and GFP for each population, median FSC and GFP for each gated population
  #Titir & Julie

  #frequency of cells inside each gate
#  cyto_stats_compute(transformed_timepoint_gating_set, #frequency of cells inside each gate
#                                      parent = c("Single_cells"),
#                                     alias = c("zero_copy", "one_copy", "two_or_more_copy"),
#                                      stat="freq",
#                                      save_as = paste0("01_02_04_v2_stats_freq_",prefix,".csv") #writes to working directory
#                                      )
  #median B2-A and FSC values under the single cell gate
#  cyto_stats_compute(transformed_timepoint_gating_set,
#                                     parent = c("Single_cells"),
#                                     alias  = c("Single_cells"),
#                                     channels = c("FSC-A", "B2-A"),
#                                     stat="median",
#                                    save_as = paste0("01_02_04_v2_stats_median_overall_", prefix,".csv"))

  #cell number
#  cyto_stats_compute(transformed_timepoint_gating_set,
#                                             parent = c("Single_cells"),
#                                             alias  = c("Single_cells"),
#                                             channels = c("FSC-A", "B2-A"),
#                                             stat="count",
#                                             save_as = paste0("01_02_04_v2_stats_cell_number_", prefix,".csv"))

  #median B2-A and FSC values of cells in each gate
#  cyto_stats_compute(transformed_timepoint_gating_set,
#                                              parent = c("Single_cells"),
#                                              alias  = c("zero_copy", "one_copy", "two_or_more_copy"),
#                                              channels = c("FSC-A", "B2-A"),
#                                              stat="median",
#                                              save_as = paste0("01_02_04_v2_stats_median_gatewise_", prefix,".csv"))

  #raw transformed B2-A and FSC values for each cell (not the median)
  timepoint_raw_list <- cyto_extract(transformed_timepoint_gating_set, parent = "Single_cells", raw = TRUE, channels = c("FSC-A", "B2-A")) #raw flow data of each single cell as a list of matrices
  map_df(timepoint_raw_list, ~as.data.frame(.x), .id="name") %>% #convert to df, put list name in new column
   mutate(name = as.factor(name)) %>% #convert `name` to factor
   left_join(experiment_details %>% #join by name column to add metadata
               mutate(generation = as.factor(unique(experiment_details$generation)))) %>%
   mutate(B2A_FSC = `B2-A`/`FSC-A`) %>% #compute normalized fluor
   write_csv(paste0("01_02_04_v2_SingleCellDistributions_",prefix,".csv"))
}

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Julie

#map(folders[-1], analyze_all_exp, my_markers, gating_template = "cytek_gating.csv")
try(map(folders[5:length(folders)],analyze_all_exp, my_markers, gating_template = "cytek_gating_01_02_04_v2.csv"))
try(map(folders[23],analyze_all_exp, my_markers, gating_template = "cytek_gating_01_02_04_v2.csv"))
#STEP 7:  Combine stats_freq.csv and stats_median.csv files into a single dataframe
#Pull in all stats_* files from directories and assemble into a single dataframe
#Author: Julie

#list.files(path = ".", pattern = "01_02_04_v2_stats_freq") %>%
#  read_csv() %>%
#  write_csv(file = "01_02_04_v2_stats_freq_all_timepoints.csv")

#list.files(path = ".", pattern = "01_02_04_v2_stats_median_overall") %>%
#  read_csv() %>%
#  write_csv(file = "01_02_04_v2_stats_median_overall_all_timepoints.csv")

#list.files(path = ".", pattern = "01_02_04_v2_stats_median_gatewise") %>%
#  read_csv() %>%
#  write_csv(file = "01_02_04_v2_stats_median_gatewise_all_timepoints.csv")

#list.files(path = ".", pattern = "01_02_04_v2_stats_cell_number") %>%
#  read_csv() %>%
#  write_csv(file = "01_02_04_v2_stats_cell_number_all_timepoints.csv")

list.files(path = ".", pattern = "01_02_04_v2_SingleCellDistributions") %>%
  read_csv() %>%
  write_csv(file = "01_02_04_v2_SingleCellDistributions_all_timepoints.csv")

#STEP 8: Plot time series & assess gates
#Determine whether =>83% of controls are in the correct gate
#Make plots
#Author: Grace & Julie

# read in frequency csv, median csvs, cell numbers csvs, single cell distributions for all timepoints
freq = read_csv("01_02_04_v2_stats_freq_all_timepoints.csv") %>% rename(Gate = Population)
medians = read_csv("01_02_04_v2_stats_median_overall_all_timepoints.csv")
medians_bygate = read_csv("01_02_04_v2_stats_median_gatewise_all_timepoints.csv")
cell_numbers = read_csv("01_02_04_v2_stats_cell_number_all_timepoints.csv")
sc_distr_alltimepoints <- read.csv("01_02_04_v2_SingleCellDistributions_all_timepoints.csv", stringsAsFactors = T) %>% mutate(generation = factor(generation, levels = unique(sc_distr_alltimepoints$generation)))

# add cell number column to freq table
freq = left_join(freq, cell_numbers) %>%
  select(-Marker)

freq %>%
  filter(Count>70000) %>%
  #filter(str_detect(Description, "control"), Gate == "zero_copy") %>%
  filter(Description == "0 copy control", Gate == "one_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>%
  ggplot(aes(Frequency)) + geom_boxplot()
  #summarize(IQR(Frequency))

#determine upper threshold for zero copy gate, aka False Negative Rate of Detecting One Copy
##Everyone should be start as One Copy (atleast) except the Zero Copy Control
freq %>%
  filter(Count>70000) %>%
  #filter(str_detect(Description, "control"), Gate == "zero_copy") %>%
  filter(Description == "1 copy control", Gate == "zero_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>%
  #ggplot(aes(Frequency)) + #geom_histogram(bins = 50) #geom_boxplot()
  summarize(IQR(Frequency))
#One Copy FN Rate is  2.78% = (median+IQR)

freq %>%
  filter(Count>70000) %>%
  filter(Description == "1 copy control", Gate == "two_or_more_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>%
  #ggplot(aes(Frequency)) + geom_boxplot()
  summarize(median(Frequency) + 1.5*IQR(Frequency)) #8.58

#Determine CNV False Negative Rate. 2 copy control falling into 1copy or 0copy Gates: median + IQR
freq %>%
  filter(Count>70000) %>%
  #filter(str_detect(Description, "control"), Gate == "zero_copy") %>%
  filter(Description == "2 copy control", Gate != "two_or_more_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>%
  #ggplot(aes(Frequency)) + geom_boxplot()
summarize(IQR(Frequency))
# CNV FN Rate is 1.33% = median + IQR

#Determine lower threshold for zero copy control gate
freq %>%
  filter(Count>70000) %>%
  filter(Description == "0 copy control", Gate == "zero_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>% #View()
  #ggplot(aes(Frequency)) + geom_boxplot()
  #summarize(min(Frequency))
  summarize(median(Frequency)-1.5*IQR(Frequency))

#Determine lower threshold for one copy control gate
freq %>%
  filter(Count>70000) %>%
  filter(Description == "1 copy control", Gate == "one_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>% #View()
  #ggplot(aes(Frequency)) + geom_boxplot()
  #summarize(min(Frequency))
  summarize(median(Frequency)-1.5*IQR(Frequency)) #87.4, the left end of the lower whisker
  #summarize(1.5*IQR(Frequency))

#Determine lower threshold for two copy control gate
freq %>%
  filter(Count>70000) %>%
  filter(Description == "2 copy control", Gate == "two_or_more_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>%
  ggplot(aes(Frequency)) + geom_boxplot()
  #summarize(min(Frequency))
  summarize(median(Frequency)-1.5*IQR(Frequency)) #92.3, the left end of the lower whisker
#summarize(1.5*IQR(Frequency))
  fails = freq %>%
  filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
# check controls are in their proper gates
  filter(str_detect(Description, "control")) %>%
  select(Description, Strain, generation, Gate, Frequency, name, Count) %>%
  mutate(flag = case_when(Strain == "DGY1" & Gate == "zero_copy" & Frequency >= 95 ~ "pass",
                          Strain == "DGY1" & Gate == "zero_copy" & Frequency < 95 ~ "fail",
                          Strain == "DGY1" & Gate == "one_copy" & Frequency >= 10 ~ "fail",
                          Strain == "DGY1" & Gate == "two_or_more_copy" & Frequency >=10 ~ "fail",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency >= 79 ~ "pass",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency < 79 ~ "fail",
                          Strain == "DGY500" & Gate == "zero_copy" & Frequency >= 10 ~ "fail",
                          Strain == "DGY500" & Gate == "two_or_more_copy" & Frequency >= 10 ~ "fail",
                          Strain == "DGY1315" & Gate == "two_or_more_copy" & Frequency >= 79 ~ "pass",
                          Strain == "DGY1315" & Gate == "two_or_more_copy" & Frequency < 79 ~ "fail",
                          Strain == "DGY1315" & Gate == "zero_copy" & Frequency >= 10 ~ "fail",
                          Strain == "DGY1315" & Gate == "one_copy" & Frequency >= 10 ~ "fail"
                          ))%>%
  filter(flag == "fail") %>%
  arrange(Description)
  View(fails)
  #fails %>% write_csv("01_02_04_v2_83_fail.csv")
  #fails %>% write_csv("01_02_04_v2_fail_calc_thres_stringent_.csv")
  #fails %>% write_csv("01_02_04_v2_79_10_fail_.csv")
# plot controls over time
freq %>%
filter(Count>70000) %>%
  filter(str_detect(Description, "control")) %>%
  select(Type, Strain, Description, generation, Gate, Frequency, Count) %>%
  anti_join(fails) %>% #exclude the contaminated controls timepoints (the failed timepoints)
  ggplot(aes(generation, Frequency, color = Gate)) +
  geom_line() +
  facet_wrap(~Description) +
  ylab("% of cells in gate") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme(text = element_text(size=12))

# plot proportion of population in each gate over time for all experimental
plot_list = list()
i=1
for(exp in unique(freq$Description)) {
  #print(plot_dist(obs))
  plot_list[[i]] = freq %>%
    filter(Count>70000) %>%
    #filter(generation != 79, generation != 116,generation != 182,generation != 252) %>%
    filter(Description==exp) %>%
    ggplot(aes(generation, Frequency, color = Gate)) +
    geom_line() +
    facet_wrap(~sample) +
    ylab("% of cells in gate") +
    theme_minimal()+
    scale_x_continuous(breaks=seq(0,250,50))+
    theme(text = element_text(size=12))
  i = i+1
}
names(plot_list) = unique(freq$Description)
plot_list$`GAP1 WT architecture` # change index to view replicates for different genetic backgrounds
plot_list$`GAP1 ARS KO`
plot_list$`GAP1 LTR KO`
plot_list$`GAP1 LTR + ARS KO`

# plot proportion of the population with a CNV over time
freq %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  #filter(generation != 79, generation != 116,generation != 182,generation != 252) %>%
  group_by(sample, generation) %>%
  mutate(prop_CNV = sum(Frequency)) %>% #View()
  select(sample, generation, Description, prop_CNV) %>%
  distinct() %>%
  ggplot(aes(generation, prop_CNV, color = sample)) +
  geom_line() +
  geom_point()+
  facet_wrap(~Description) +
  ylab("Proportion of the population with GAP1 CNV") +
  scale_color_manual(values = c("#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59", #WT,5,gold
"#637EE7","#6F88E9","#7B92EA","#4463E2","#3053DF","#2246D7","#1E3FC3","#5766E6", #ALL,8,bluepurple
"#DE54B9","#E160BE","#E36CC3","#E578C8","#E885CD","#DB41B2","#D72FAA", #ARS,7,pink
"#54DE79","#41DB6A","#2FD75C","#26CA52","#23B84B","#60E182","#6CE38C","#78E595" #LTR,8,green
   )) +
  #scale_color_manual(values = c()
  theme_minimal() +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme(text = element_text(size=12),legend.position = "none")

#Plot proportion of the populations with a CNV over time (collapse the replicates)
freq %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  #ggplot(aes(Frequency)) + geom_histogram(bins=50) #right skewed distribution
  ggplot(aes(generation, Frequency, color = Description)) +
  scale_color_manual(values=c("#DE54B9", "#5474DE", "#54DE79", "#DEBD52"))+
  scale_fill_manual(values=c("#DE54B9", "#5474DE", "#54DE79", "#DEBD52"))+
  #geom_point(alpha = 0.5, size =1) +
  stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill=Description), alpha=0.3) +
  facet_wrap(~Description) +
  theme_minimal() +
  ggtitle("loess regression of the populations") +
  ylab("Proportion of the population with GAP1 CNV") +
  scale_x_continuous(breaks=seq(0,250,50)) +
  scale_y_continuous(breaks=seq(0,100,25))+f
  theme(text = element_text(size=14), legend.position = "none")

#plot ridgeplots of controls over time
  #sc_distr_alltimepoints %>%
  #  mutate(generation = factor(generation, levels = unique(sc_distr_alltimepoints$generation)))
  #zero <- sc_distr_alltimepoints %>% filter(Description == "0 copy control") %>%
  #write_csv(file = "sc_distributions_0copyControl_all_timepoints.csv")
  zero = read.csv("sc_distributions_0copyControl_all_timepoints.csv", stringsAsFactors = T)
  zero = zero %>% mutate(generation = factor(generation, levels = unique(zero$generation))) #convert generation to factor
  #one <- sc_distr_alltimepoints %>% filter(Description == "1 copy control") %>%
  #write_csv(file = "sc_distributions_1copyControl_all_timepoints.csv")
  one = read.csv("sc_distributions_1copyControl_all_timepoints.csv", stringsAsFactors = T)
  one = one %>% mutate(generation = factor(generation, levels = unique(one$generation)))
  #two <- sc_distr_alltimepoints %>% filter(Description == "2 copy control") %>%
  # write_csv(file = "sc_distributions_2copyControl_all_timepoints.csv")
  two = read.csv("sc_distributions_2copyControl_all_timepoints.csv", stringsAsFactors = T)
  two = two %>% mutate(generation = factor(generation, levels = unique(two$generation)))
zero_ridges = ggplot(zero, aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
    geom_density_ridges_gradient(scale = 1.0, rel_min_height = 0.01) +
    xlab("normalized fluorescence") +
    ylab("generation") +
    ggtitle("zero copy control") +
    theme_classic() +
    #scale_x_continuous("Normalized Fluorescence", limits=c(0.05,1.0), expand = c(0.01, 0), breaks = c(0.1, 0.55, 1.0)) +
    scale_y_discrete(expand = expansion(add = c(0.2, 1.0))) + #expands the graph space or else the top is cut off
    scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
    theme(
      legend.text = element_text(family="Arial", size = 12),#edit legend text font and size
      legend.title = element_blank(), #remove legend title
      legend.position = 'none', #remove the legend
      axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
      axis.text.y = element_text(family="Arial", size = 10, color = "black")
    )
one_ridges = ggplot(one, aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
  geom_density_ridges_gradient(scale = 1.0, rel_min_height = 0.01) +
  #xlab("normalized fluorescence") +
  ylab("generation") +
  ggtitle("one copy control") +
  theme_classic() +
  scale_x_continuous(limits=c(0.0,2.5), breaks = c(0, 1, 2, 2.5)) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.0))) + #expands the graph space or else the top is cut off
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
  theme(
    legend.position = 'none', #remove the legend
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
one_ridges

two_ridges = ggplot(two, aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
  geom_density_ridges_gradient(scale = 1.0, rel_min_height = 0.01) +
  #xlab("normalized fluorescence") +
  ylab("generation") +
  ggtitle("two copy control") +
  theme_classic() +
  scale_x_continuous(limits=c(0.0,2.5), breaks = c(0, 1, 2, 2.5)) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.0))) + #expands the graph space or else the top is cut off
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
  theme(
    legend.position = 'none', #remove the legend
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
two_ridges

cowplot::plot_grid(zero_ridges, one_ridges, two_ridges,labels = c("A"), nrow=1, align = "h")

controls = bind_rows(zero, one, two)

#all Controls in one ggplot, and facet by Description
ggplot(controls, aes(B2A_FSC, generation, fill = Description)) +
  geom_density_ridges(scale = 1) +
  facet_grid(~Description) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.0)))+
  #scale_fill_brewer(type = "seq", palette = 5, direction = 1) +
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(4, "Greens")[-1])) +
  scale_x_continuous("normalized fluorescence", limits=c(0, 2.5), breaks = c(0, 1, 2, 2.5), labels = c(0,1,2,2.5)) +
  theme_classic() +
  theme(
      legend.position = 'none', #remove the legend
      axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
      axis.text.y = element_text(family="Arial", size = 10, color = "black"),
      strip.background = element_blank(), #removed box around facet title
      strip.text = element_text(size=12)
  )
ggsave("controls_generation_ridgeplot_Facet.png")

##############################################
# STEP 9:  Quantify CNV dynamics (Lauer et al. 2018)
# Author: Julie
  # 1) First, calculate Tup, the generation at which CNVs are initially detected, (Lang et al. 2011 and Lauer et al. 2018)
  # To do that, calculate the false positive rate for CNV detection (threshold), which I will define as the median frequency of 1 copy control cells appearing in the two_copy_or_more gate and gate across generations 8-260 plus the interquartile range (IQR) if the distribution of one copy controls appearing in the CNV gate as NOT normal. If the distribution is normal, then use the mean plus one standard deviation like in Lauer et al. 2018. Like in Lauer et al. 2018, samples surpassing this threshold is considered to contain CNVs.

#CNV False Positive Rate is defined by the frequency of the 1 copy control strain appearing in the CNV gate which is the Two_or_more copy gate.
CNV_false_pos_df = freq %>%
  filter(Count>70000) %>%
  anti_join(fails) %>%
  filter(Type == "1_copy_ctrl") %>%
  filter(Gate %in% c("two_or_more_copy")) %>%
  select(Type, Strain, Description, generation, Gate, Frequency, Count)

# draw a histogram, see distribution is normal by eye
hist(CNV_false_pos_df$Frequency) #looks normal but could be left skewed
  abline(v = median(CNV_false_pos_df$Frequency),col = "red",lwd = 1.5)
  abline(v = mean(CNV_false_pos_df$Frequency), col = "blue", lwd = 1.5)

#Test for normality
shapiro.test(CNV_false_pos_df$Frequency) #null hypothesis is that the distribution is normal. if p <0.05, then it rejects the null hypothesis and so the distribution is NOT normal.
    #W = 0.95715, p-value = 0.4886
# Our distribution is normal. Use the mean + 1SD as the threshold value.
thres_mean = mean(CNV_false_pos_df$Frequency) + sd(CNV_false_pos_df$Frequency) #6.615341

#Determine Tup for each population, the generation when CNVs first appear.
Tup_per_pop = freq %>%
  filter(Count>70000) %>%
  anti_join(fails) %>%
  filter(Type == "Experimental", Gate == "two_or_more_copy", Frequency >= thres_mean) %>%
  select(Type, Strain, Description, sample, generation, Gate, Frequency) %>% #View()
  group_by(sample) %>%
  slice(which.min(generation))

ggplot(Tup_per_pop, aes(reorder(Description, -generation),generation, fill = Description)) +
  geom_boxplot() +
  xlab("genotype") +
  #scale_fill_manual(values=c("#00B0F6", "#00BF7D", "#E76BF3", "#B79F00")) +  #c(ARS,LTR+ARS,LTR,WT)
  scale_fill_manual(values=c("#DE54B9", "#5474DE", "#54DE79", "#DEBD52")) +  #c(ARS,LTR+ARS,LTR,WT)
  #scale_fill_manual(values=c("#AA2787", "#2745AA", "#26AA49", "#AA8C27")) +  #c(ARS,LTR+ARS,LTR,WT)
  #scale_fill_manual(values=c(awtools::spalette))+
  ylab("generation of first CNV appearance") +
  scale_y_continuous(breaks=c(8,20, 30, 40, 50, 60, max(Tup_per_pop$generation)))+
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=12))

# ANOVA to test for significance
# One Way ANOVA because there is only 1 independent variable, genotype.
# Anova Tut: https://www.scribbr.com/statistics/anova-in-r/
Tup_anova = aov(generation~Description, data = Tup_per_pop)
summary(Tup_anova)
#              Df Sum Sq Mean Sq F value Pr(>F)
# Description  3   2377   792.3   2.569  0.078 .
# Residuals   24   7402   308.4
# Conclusion: There is no difference in the means. Genotype has no significant impact on Tup.

########### GFP over time plots
#Using outputted `medians` .csv files, plot median GFP fluorescence normalized over median FSC-A over time for experimental
  #overlay 0,1,2 controls on same graph as experimental with gray lines
  norm_medians = medians %>%
    pivot_wider(names_from = Marker, values_from = MedFI) %>%
    mutate(NormMedGFP = GFP/`FSC-A`) #normalization

  #Mutate description of controls so we can graph them on experimental facet plots
  dups = norm_medians %>% arrange(Description) %>%
  slice(rep(1:56, each = 4)) %>% #repeat each control row 4 times because was have 4 facetplots
  mutate(Description = rep(c("GAP1 ARS KO", "GAP1 LTR + ARS KO", "GAP1 LTR KO","GAP1 WT architecture"), times=56))

  #merge back to experimental rows from norm_medians df
  adj_norm_medians = merge(norm_medians %>% filter(Type == "Experimental"), dups, all = TRUE) %>% #merge back to experimental rows
  arrange(Description, generation)

#Clean up data frame. Remove timepoints/controls that are abnormal.
  #I don't think removal of timepoints is justified. Our data is very variable, so nothing is an outlier if you draw a boxplot of all GFP or Freq from ONE generation from ONE Genotype. Spikes that we see in one timepoint could be real since we sample every 8-10 generations. Lineages containing CNVs can indeed increase in that time frame, especially since we see CNVs starting at gen8.
clean_adj_norm_medians = adj_norm_medians %>%
  anti_join(adj_norm_medians %>% filter(generation == 182 & Type == "1_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 203 & Type == "1_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 231 & Type == "0_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 260 & Type == "1_copy_ctrl"))
#anti_join() gap1_4 two copy gate samples at g79, 124, 231, 252

#Graph experimental with along controls
  ggplot(clean_adj_norm_medians, aes(generation, NormMedGFP, color= sample)) +
  geom_line(aes(linetype = Type)) +
  scale_linetype_manual(values = c("dashed", "dashed", "dashed", "solid")) +
  scale_color_manual(values=c("gray", "gray", "gray", #controls
                              rep("orange", 5),  #WT,5
                              rep("purple", 8),  #ALL,8
                              rep("blue", 7), #ARS, 7
                              rep("pink", 8)  #LTR,8
                              ))+
  facet_wrap(~Description) +
  ylab("normalized median fluorescence") +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme_classic()
#  + theme(legend.position = "none",
#        text = element_text(size=12))

# Combine the replicates/populations and plot median normalized fluorescence over time
# dashed gray controls lines on top of the experiment lineplot
ggplot(data = clean_adj_norm_medians %>% filter(Type == "Experimental"),
      mapping = aes(generation, NormMedGFP, color = Description)) +
      stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill = Description), alpha=0.3) + #experimentals loess regression with standard error cloud
    geom_line(mapping = aes(generation, NormMedGFP, color = Type, linetype = Type),
                data = clean_adj_norm_medians %>% filter(Type != "Experimental")) + #controls lineplot
    scale_linetype_manual(values = c("dashed", "dashed", "dashed"))+ #dashed lines for controls
    scale_color_manual(values=c("gray", "gray", "gray","#DE54B9", "#5474DE", "#54DE79", "#DEBD52"))+ #line color
    scale_fill_manual(values=c("#DE54B9", "#5474DE", "#54DE79", "#DEBD52"))+ #standard error cloud color
    facet_wrap(~Description) +
    theme_classic() +
    ggtitle("loess regression of combined populations") + theme(plot.title = element_text(hjust = 0.5)) +
    ylab("normalized median fluorescence") +
    scale_x_continuous(breaks=seq(0,250,50)) +
    theme(text = element_text(size=12),legend.position = "none") #I did it!

  #raw - standalone.
  adj_norm_medians %>%
    filter(Type != "Experimental") %>%
    ggplot(aes(generation, NormMedGFP, color = Type))+
    geom_line(aes(linetype = Type)) +
    scale_linetype_manual(values = c("dashed", "dashed", "dashed"))+
    scale_color_manual(values = c("gray", "gray", "gray"))

  #geom_line() without ggplot()
  geom_line(mapping = aes(generation, NormMedGFP, color = Type, linetype = Type),
            data = adj_norm_medians %>% filter(Type != "Experimental"))+
    scale_linetype_manual(values = c("dashed", "dashed", "dashed"))+
    scale_color_manual(values = c("gray", "gray", "gray"))

#g231 0copy control looks abnormal. Investigate and see if whole timepoint should be thrown out
  # GFP values are abnormally high, FSC-A values are abnormally low.
  # I can do some hypothesis testing to see they are not part of the null GFP distribution or null FSC distribution.
  adj_norm_medians %>% filter(Type == "0_copy_ctrl", generation >= 200) %>% View()

#g203 1copy control look abnormal. GFP values are abnormally low. FSC-A values are normal.
  adj_norm_medians %>% filter(Type == "1_copy_ctrl", generation >= 160 & generation <= 211) %>% View()

#MY PALLETEE
#GOLD - #DEBD52
#Soft Blue - #54A2DE
#Soft Blue Purple - #5474DE
#Soft Purple - #6254DE
#Hot Turquoise - #54DEBE
# Hot Pink - #BE54DE
# Light Green - #74DE54

#Square
#DEBD52 #gold
#54DE79 #neon green
#5474DE #purple or #Soft Blue Purple - #5474DE
#DE54B9 #pink

#Square
#26AA49 #true green
#2745AA #purple
#AA2787 #magent
#AA8C27 #golde


#Gold Metallic (6)
#DEBD52
#DBB741
#D7B02F
#CAA426
#D9BB59
#B89523

#Green Machalite (8)
#54DE79
#41DB6A
#2FD75C
#26CA52
#23B84B
#60E182
#6CE38C
#78E595

#royal blue light (8)
#5774E5
#637EE7
#6F88E9
#7B92EA
#4463E2
#3053DF
#2246D7
#1E3FC3
#5766E6

#Pink (8)
#DE54B9
#E160BE
#E36CC3
#E578C8
#E885CD
#DB41B2
#D72FAA
#BE54DE #medium orchid
#B741DB
#B02FD7
#A426CA

  #ARS = BLUE  #00B0F6
  # LTR = PINK  #E76BF3
  # GREEN #00BF7D or #7CAE00
  # GOLD = #B79F00 or #C49A00

