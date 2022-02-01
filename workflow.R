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

#use flowWorkspace::pData to annotate the experiment details file associated with the gating set
experiment_details <- read_csv(exp_details_path) #import experiment-details.csv
for(i in 1:length(names(experiment_details))){
  flowWorkspace::pData(timepoint_gating_set)[names(experiment_details[i])]<-experiment_details[i]
  }

#file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv") #rename the experiment-markers.csv file. Need to do once.

#STEP 3:  Perform gating on gating set
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.
#Author: Titir & Julie

#transform the data
# looks useful if I want to choose different transformation: https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR-Transformations.html
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

# note:if you already have a gating template and don't need to draw gates, then skip cyto_draw, use cyto_gatingTemplate_apply to apply the gating template.csv to your gating set
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

#STEP 4:  Generate statistics tables and single cell data tables
#Author: Titir & Julie

stats_freq_01 <- cyto_stats_compute(transformed_timepoint_gating_set01,
                                  parent = "Single_cells",
                                  alias = c("zero_copy", "one_copy", "two_or_more_copy"),
                                  stat="freq",
                                  save_as = "stats_freq_01.csv")

timepoint_raw_list <- cyto_extract(transformed_timepoint_gating_set, parent = "Single_cells", raw = TRUE, channels = c("FSC-A", "B2-A")) #raw flow data of each single cell as a list of matrices

map_df(timepoint_raw_list, ~as.data.frame(.x), .id="name") %>% #convert to df, put list name in new column
  mutate(name = as.factor(name)) %>% #convert `name` to factor
  left_join(experiment_details %>% #join by name column to add metadata
  mutate(generation = as.factor(unique(experiment_details$generation)))) %>%
  mutate(B2A_FSC = `B2-A`/`FSC-A`) %>% #compute normalized fluor
  write_csv(paste0("01_02_04_v2_SingleCellDistributions_",prefix,".csv"))

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
experiment_details
  #3. specify markers for that gating set
  markernames(timepoint_gating_set)<-my_markers

  #4. transform data
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

  #5. apply gating-template.csv to transformed gating set
  cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= gating_template)
#  cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= "cytek_gating_01_02_04_v2.csv")

  #6. write stats: freq file for % of cells inside each gate, median FSC and GFP for each population, median FSC and GFP for each gated population
  #Titir & Julie

  #frequency of cells inside each gate
#  cyto_stats_compute(transformed_timepoint_gating_set, #frequency of cells inside each gate
#                                      parent = c("Single_cells"),
#                                     alias = c("zero_copy", "one_copy", "two_or_more_copy"),
#                                      stat="freq",
#                                      save_as = paste0("01_02_04_v2_stats_freq_",prefix,".csv") #writes to working directory
#                                      )

  #cell number
#  cyto_stats_compute(transformed_timepoint_gating_set,
#                                             parent = c("Single_cells"),
#                                             alias  = c("Single_cells"),
#                                             stat="count",
#                                             save_as = paste0("01_02_04_v2_stats_cell_number_", prefix,".csv"))
  #get cell count from each gate
  #cyto_stats_compute(transformed_timepoint_gating_set,
  #                   parent = "Single_cells",
  #                  alias = c("Single_cells","zero_copy", "one_copy", "two_or_more_copy"),
  #                   stat="count"
  #                   save_as = paste0("01_02_04_v2_stats_cell_number_in_Gates", prefix,".csv"))
  ##cyto_stat_compute() computes freq and count values that DO NOT match % in cyto_plot_gating_scheme. Instead, use gs_pop_get_stats()
  #get cell count from each gate
  gs_pop_get_stats(transformed_timepoint_gating_set, c("Single_cells", "zero_copy", "one_copy", "two_or_more_copy")) %>%
    rename(Gate = pop, name = sample, Count = count) %>%
    left_join(experiment_details) %>%
    write_csv(paste0("01_02_04_v2_fw_counts_", prefix,".csv"))
  #get frequency of cells inside each gate
  gs_pop_get_stats(transformed_timepoint_gating_set, c("Single_cells","zero_copy", "one_copy", "two_or_more_copy"), type = "percent") %>%
    rename(Gate = pop, name = sample, Frequency = percent) %>%
    left_join(experiment_details) %>%
    write_csv(paste0("01_02_04_v2_fw_freq_", prefix,".csv"))

  #raw transformed B2-A and FSC values for each cell (not the median)
#  timepoint_raw_list <- cyto_extract(transformed_timepoint_gating_set, parent = "Single_cells", raw = TRUE, channels = c("FSC-A", "B2-A")) #raw flow data of each single cell as a list of matrices
#  map_df(timepoint_raw_list, ~as.data.frame(.x), .id="name") %>% #convert to df, put list name in new column
#   mutate(name = as.factor(name)) %>% #convert `name` to factor
#   left_join(experiment_details %>% #join by name column to add metadata
#               mutate(generation = as.factor(unique(experiment_details$generation)))) %>%
#   mutate(B2A_FSC = `B2-A`/`FSC-A`) %>% #compute normalized fluor for each cell
#   write_csv(paste0("01_02_04_v2_SingleCellDistributions_",prefix,".csv"))
}

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Julie

#map(folders[-1], analyze_all_exp, my_markers, gating_template = "cytek_gating.csv")
try(map(folders[4:length(folders)],analyze_all_exp, my_markers, gating_template = "cytek_gating_01_02_04_v2.csv"))
try(map(folders[23],analyze_all_exp, my_markers, gating_template = "cytek_gating_01_02_04_v2.csv"))
#STEP 7:  Combine stats_freq.csv files into a single dataframe
#Pull in all stats_* files from directories and assemble into a single dataframe
#Author: Julie

#list.files(path = ".", pattern = "01_02_04_v2_stats_freq") %>%
#  read_csv() %>%
#  write_csv(file = "01_02_04_v2_stats_freq_all_timepoints.csv")

#list.files(path = ".", pattern = "01_02_04_v2_stats_cell_number") %>%
#  read_csv() %>%
#  write_csv(file = "01_02_04_v2_stats_cell_number_all_timepoints.csv")

list.files(path = ".", pattern = "01_02_04_v2_fw_counts_([0-9])+_EE_GAP1_ArchMuts_2021") %>%
  read_csv() %>%
  write_csv(file = "01_02_04_v2_fw_counts_all_timepoints.csv")

list.files(path = ".", pattern = "01_02_04_v2_fw_freq_([0-9])+_EE_GAP1_ArchMuts_2021") %>%
  read_csv() %>%
  write_csv(file = "01_02_04_v2_fw_freq_all_timepoints.csv")

#do on hpc because large files
#list.files(path = ".", pattern = "01_02_04_v2_SingleCellDistributions") %>%
#  read_csv() %>%
#  write_csv(file = "01_02_04_v2_SingleCellDistributions_all_timepoints.csv")

#STEP 8: Plot ridgeplots, time series, & assess gates
#Determine whether =>83% of controls are in the correct gate
#Make plots
#Author: Grace & Julie

# read in frequency csv, cell numbers csvs, single cell distributions for all timepoints
freq = read_csv("01_02_04_v2_stats_freq_all_timepoints.csv") %>% rename(Gate = Population)
fw_freq = read_csv("01_02_04_v2_fw_freq_all_timepoints.csv") #%>% rename(Frequency = frequency)
cell_numbers = read_csv("01_02_04_v2_stats_cell_number_all_timepoints.csv")
fw_counts= read_csv("01_02_04_v2_fw_counts_all_timepoints.csv")
sc_distr_alltimepoints <- read.csv("01_02_04_v2_SingleCellDistributions_all_timepoints.csv", stringsAsFactors = T) %>% mutate(generation = factor(generation, levels = unique(generation)))
freq = left_join(freq, cell_numbers) %>% # add cell number column to freq table
  select(-Marker)

fw_freq_and_counts =
  fw_counts %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(fw_freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate)

#determine upper threshold for zero copy gate, aka False Negative Rate of Detecting One Copy
##Everyone should be start as One Copy (atleast) except the Zero Copy Control
# one copy control falling into the zero copy or two+ copy gates
freq %>%
  filter(Count>70000) %>%
  #filter(str_detect(Description, "control"), Gate == "zero_copy") %>%
  filter(Description == "1 copy control" & Gate == "zero_copy" | Description == "1 copy control" & Gate == "two_or_more_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>% #View()
  #ggplot(aes(Frequency)) +
  #geom_histogram(bins = 50)
  #  geom_boxplot()
  summarize(IQR(Frequency)+median(Frequency)) #One Copy FN Rate is  6.86% = (median+IQR)

freq %>%
  filter(Count>70000) %>%
  filter(Description == "1 copy control", Gate == "two_or_more_copy") %>%
  select(Type, Strain, Description, sample,generation, Gate, Frequency, Count) %>%
  #ggplot(aes(Frequency)) + geom_boxplot()
  #summarize(median(Frequency) + 1.5*IQR(Frequency)) #8.58
  summarize(median(Frequency) + IQR(Frequency))

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

#Table of low cell observations, convenient to have to anti_join() in further steps
freq %>% filter(Count <7000) %>% View()

lowcell = freq %>%
  filter(Count <7000) %>%
  mutate(generation = factor(generation, levels = unique(generation))) %>%
  select(-Count)

lowcell = fw_freq_and_counts %>%
  filter(Count <7000) %>%
  mutate(generation = factor(generation, levels = unique(generation))) %>% #View()
  select(-Count)

#df %>% anti_join(lowcell) #example code to remove lowcell df from bigger table

## check controls are in their proper gates
  fails = fw_freq_and_counts %>%
    #fails = freq %>%
  filter(Count>70000) %>% # exclude any well/timepoint with less than 70,000 single cells
  filter(str_detect(Description, "control")) %>%
  select(Description, Strain, generation, Gate, Frequency, name, Count) %>%
  mutate(flag = case_when(Strain == "DGY1" & Gate == "zero_copy" & Frequency >= 95 ~ "pass",
                          Strain == "DGY1" & Gate == "zero_copy" & Frequency < 95 ~ "fail",
                          Strain == "DGY1" & Gate == "one_copy" & Frequency >= 10 ~ "fail",
                          Strain == "DGY1" & Gate == "two_or_more_copy" & Frequency >=11 ~ "fail",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency >= 79 ~ "pass",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency < 79 ~ "fail",
                          Strain == "DGY500" & Gate == "zero_copy" & Frequency >= 11 ~ "fail",
                          Strain == "DGY500" & Gate == "two_or_more_copy" & Frequency >= 11 ~ "fail",
                          Strain == "DGY1315" & Gate == "two_or_more_copy" & Frequency >= 79 ~ "pass",
                          Strain == "DGY1315" & Gate == "two_or_more_copy" & Frequency < 79 ~ "fail",
                          Strain == "DGY1315" & Gate == "zero_copy" & Frequency >= 11 ~ "fail",
                          Strain == "DGY1315" & Gate == "one_copy" & Frequency >= 11 ~ "fail"
                          ))%>%
  filter(flag == "fail") %>%
  arrange(Description)
  View(fails)
  #fails %>% write_csv("01_02_04_v2_83_fail.csv")
  #fails %>% write_csv("01_02_04_v2_fail_calc_thres_stringent_.csv")
  #fails %>% write_csv("01_02_04_v2_79_10_fail_.csv")
  fails %>% write_csv("01_02_04_v2_fw_79_11_fail.csv")

# plot controls over time
freq %>%
filter(Count>70000) %>%
  filter(str_detect(Description, "control")) %>%
  select(Type, Strain, Description, generation, Gate, Frequency, Count) %>%
  #anti_join(fails) %>% #exclude the contaminated controls timepoints (the failed timepoints)
  ggplot(aes(generation, Frequency, color = Gate)) +
  geom_line() +
  facet_wrap(~Description) +
  ylab("% of cells in gate") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme(text = element_text(size=12))

# plot controls over time
fw_freq_and_counts %>%
  filter(Count>70000,
          str_detect(Description, "control")) %>%
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

#dev
plot_list = list()
i=1
for(exp in unique(fw_freq_and_counts$Description)) {
  #print(plot_dist(obs))
  plot_list[[i]] = fw_freq_and_counts %>%
    #filter(Gate == "Single_cells" & Count>70000) %>%
    #filter(!(Gate == "Single_cells"))%>%
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
names(plot_list) = unique(fw_freq_and_counts$Description)
plot_list$`GAP1 WT architecture` # change index to view replicates for different genetic backgrounds
plot_list$`GAP1 ARS KO`
plot_list$`GAP1 LTR KO`
plot_list$`GAP1 LTR + ARS KO`

#Plot ridgeplots of each population to catch contaminated timepoints/outlier values.
  #(Later will calculate Sup values for each population)
#I could write a Function that plots ridgeplots. Then I can write use map() to apply this FUNCtion to all population.csv since I have 32 pops
pop_files = list.files(pattern = "sc_distributions_")[-1:-3]
#file_name = pop_files[9] #dev
make_ridgeplots = function(file_name){
  pop_name = sub("sc_distributions_", "", sub("_all_timepoints.csv","", file_name))

  pop_data = read.csv(file_name, stringsAsFactors = T) %>%
    mutate(generation = factor(generation, levels = unique(generation)))

  lowcell = count(pop_data, generation, sample) %>% filter(n < 70000)

  pop_data %>%
    anti_join(lowcell) %>%
    ggplot(aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
    geom_density_ridges_gradient(scale = 2.0, rel_min_height = 0.01) +
    xlab("Normalized fluorescence (a.u.)") +
    ylab("Generation") +
    ggtitle(paste0(pop_name)) +
    theme_classic() +
    scale_x_continuous(limits=c(0.0,3), breaks = c(0, 1, 2, 3.0)) +
    scale_y_discrete(expand = expansion(add = c(0.2, 2.5))) + #expands the graph space or else the top is cut off
    scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
    theme(
      legend.position = 'none', #remove the legend
      axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
      axis.text.y = element_text(family="Arial", size = 10, color = "black")
    )
  ggsave(paste0(pop_name,"_ridgeplot_scale2.png"))
}

map(pop_files, make_ridgeplots)



# plot proportion of the population with a CNV over time
#freq %>%
fw_freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(fails)  %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  group_by(sample, generation) %>%
  mutate(prop_CNV = sum(Frequency)) %>% #View()
  select(sample, generation, Description, prop_CNV) %>%
  distinct() %>%
  ggplot(aes(generation, prop_CNV, color = sample)) +
  geom_line() +
  #geom_point()+
  facet_wrap(~Description) +
  xlab("Generation") +
  ylab("Proportion of the population with GAP1 CNV") +
  scale_color_manual(values = c("#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59", #WT,5,gold
"#637EE7","#6F88E9","#7B92EA","#4463E2","#3053DF","#2246D7","#1E3FC3","#5766E6", #ALL,8,bluepurple
"#DE54B9","#E160BE","#E36CC3","#E578C8","#E885CD","#DB41B2","#D72FAA", #ARS,7,pink
"#54DE79","#41DB6A","#2FD75C","#26CA52","#23B84B","#60E182","#6CE38C","#78E595" #LTR,8,green
   )) +
  #scale_color_manual(values = c()
  theme_classic() +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme(text = element_text(size=12),
        legend.position = "none",
        axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
        axis.text.y = element_text(family="Arial", size = 10, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=12)
        )

#Plot proportion of the populations with a CNV over time (collapse the replicates)
#freq %>%
fw_freq_and_counts %>%
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
  scale_y_continuous(breaks=seq(0,100,25))+
  theme(text = element_text(size=12),
        legend.position = "none",
        axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
        axis.text.y = element_text(family="Arial", size = 10, color = "black"))


###################################################
# plot ridgeplots (histograms):
# instead of looking at the MEDIAN GFP values per population per generation,
# we want to see the ENTIRE DISTRIBUTION of GFP of the cells per population per generation.

#plot ridgeplots of controls over time
#sc_distr_alltimepoints %>%
#mutate(generation = factor(generation, levels = unique(generation))) %>%
#filter(Description == "0 copy control") %>%
#write_csv(file = "sc_distributions_0copyControl_all_timepoints.csv")
zero = read.csv("sc_distributions_0copyControl_all_timepoints.csv", stringsAsFactors = T) %>%
  mutate(generation = factor(generation, levels = unique(generation))) #convert generation to factor
#one <- sc_distr_alltimepoints %>% filter(Description == "1 copy control") %>%
  #write_csv(file = "sc_distributions_1copyControl_all_timepoints.csv") #do once
one = read.csv("sc_distributions_1copyControl_all_timepoints.csv", stringsAsFactors = T) %>%
    mutate(generation = factor(generation, levels = unique(generation)))
#two <- sc_distr_alltimepoints %>% filter(Description == "2 copy control") %>%
  # write_csv(file = "sc_distributions_2copyControl_all_timepoints.csv") #do once
two = read.csv("sc_distributions_2copyControl_all_timepoints.csv", stringsAsFactors = T) %>%
    mutate(generation = factor(generation, levels = unique(generation)))

#all Controls in one ggplot, and facet by Description
controls = bind_rows(zero, one, two) #do once
controls %>%
  filter(!(Description == "1 copy control" & generation == 203 |
         Description == "0 copy control" & generation == 231)) %>%
ggplot(aes(B2A_FSC, generation, fill = Description)) +
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
ggsave("controls_generation_ridgeplot_excludeLowCellSamples.png")
ggsave("controls_generation_ridgeplot_Facet.png")


###### Plot normalized median mCitrine fluorescence over time
#overlay 0,1,2 controls on same graph as experimental with gray lines
###### 01-05-22 It's not accurate to normalize after taking the median GFP and median FSC-A values.
# In fact the graphs would look different. Instead use the single cell data, normalize B2-A by FSC first, then take the median of the normalized fluorescence for these line plots.
# Lauer et al. 2018 did this, see Methods - Flow Cytometry Sampling & Analysis
# For each unique `sample` calculate the median B2-A/FSC-A at each generation.

# on hpc, do once
sc_distr_alltimepoints %>%
  group_by(sample, generation) %>%
  mutate(Med_B2A_FSC = median(B2A_FSC)) %>%
  distinct(Med_B2A_FSC, .keep_all = T) %>%
  select(-FSC.A, -B2.A, -B2A_FSC) %>%
  write_csv("medians_normalized_fluor_alltimepoints.csv")

norm_medians = read_csv("medians_normalized_fluor_alltimepoints.csv") %>%
  left_join(cell_numbers) %>%
  select(-Marker) %>%
  rename(Gate = Population)

#Rename description of controls so we can graph them on experimental facet plots
relabel_controls = norm_medians %>% arrange(Description) %>%
  slice(rep(1:sum(str_detect(norm_medians$Type, "ctrl")), each = 4)) %>% #repeat each control row 4 times because was have 4 facetplots
  mutate(Description = rep(c("GAP1 ARS KO", "GAP1 LTR + ARS KO", "GAP1 LTR KO","GAP1 WT architecture"), times=sum(str_detect(norm_medians$Type, "ctrl"))))

#merge back to experimental rows from norm_medians df
adj_norm_medians = merge(norm_medians %>% filter(Type == "Experimental"), relabel_controls, all = TRUE) %>% #merge back to experimental rows
  arrange(Description, generation) %>%
  filter(Count>70000) #exclude observations with <70,000 cells

#Clean up data frame. Remove timepoints of controls that are abnormal as informed by ridgeplots
clean_adj_norm_medians = adj_norm_medians %>%
  #remove select controls timepoints based on ridgeplots
  anti_join(adj_norm_medians %>% filter(generation == 116 & Type == "2_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 182 & Type == "1_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 203 & Type == "1_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 231 & Type == "0_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 260 & Type == "1_copy_ctrl"))
#anti_join() gap1_4 two copy gate samples at g79, 124, 231, 252 -- why? not justified.

#Graph experimental with along controls
clean_adj_norm_medians %>%
  filter(!(Med_B2A_FSC<1.5 & Type == "Experimental")) %>%#filter out outliers (likely resulting from contamination) as defined by Fluor <1.5
ggplot(aes(generation, Med_B2A_FSC, color= sample)) +
  geom_line(aes(linetype = Type)) +
  scale_linetype_manual(values = c("dashed", "dashed", "dashed", "solid")) +
  scale_color_manual(values=c("gray", "gray", "gray", #controls
                              "#DEBD52", "#DBB741", "#D7B02F", "#CAA426","#D9BB59", #WT,5, golds
                              #rep("#5474DE", 8),  #ALL,8, blue/purple "#5474DE"
                              "#5774E5","#637EE7", "#6F88E9","#7B92EA","#4463E2","#3053DF","#2246D7","#1E3FC3", #ALL,8, blue/purple "#5474DE"
                              "#DE54B9","#E160BE","#E36CC3","#E578C8","#E885CD","#DB41B2","#D72FAA", #rep("#DE54B9", 7), #ARS, 7,  pink "#DE54B9"
                              "#54DE79","#41DB6A","#2FD75C","#26CA52","#23B84B","#60E182","#6CE38C","#78E595" #rep("#54DE79", 8)  #LTR,8, green "#54DE79"
  ))+
  facet_wrap(~Description) +
  xlab("Generation") +
  ylab("Median normalized fluorescence (a.u.)") +
  scale_x_continuous(breaks=seq(0,260,50)) +
  #ylim(c(1.5,2.5))+
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=12),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=12),
        axis.text.x = element_text(family="Arial", size = 12, color = "black"), #edit x-tick labels
        axis.text.y = element_text(family="Arial", size = 12, color = "black"))
ggsave("MedNormFluo_FacetPlots_NoOutliers_010722.png")
ggsave("MedNormFluro_v2_010722.png")

# Combine the replicates/populations and plot median normalized fluorescence over time
      # dashed gray controls lines on top of the experiment lineplot
clean_adj_norm_medians %>%
  filter(Type == "Experimental") %>%
  filter(!(Med_B2A_FSC<1.5 & Type == "Experimental")) %>% #filter out outliers (likely resulting from contamination)
ggplot(mapping = aes(generation, Med_B2A_FSC, color = Description)) +
  stat_smooth(method="loess", span=0.1, se=TRUE, aes(fill = Description), alpha=0.3) + #experimentals loess regression with standard error cloud
  geom_line(mapping = aes(generation, Med_B2A_FSC, color = Type, linetype = Type),
            data = clean_adj_norm_medians %>% filter(Type != "Experimental")) + #controls lineplot
  scale_linetype_manual(values = c("dashed", "dashed", "dashed"))+ #dashed lines for controls
  scale_color_manual(values=c("gray", "gray", "gray","#DE54B9", "#5474DE", "#54DE79", "#DEBD52"))+ #line color
  scale_fill_manual(values=c("#DE54B9", "#5474DE", "#54DE79", "#DEBD52"))+ #standard error cloud color
  facet_wrap(~Description) +
  theme_classic() +
  ggtitle("loess regression of combined populations") + theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Generation") +
  ylab("Median normalized fluorescence (a.u.)") +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme(legend.position = "none",
        text = element_text(size=12),
        strip.background = element_blank(), #remove box around facet title
        strip.text = element_text(size=12),
        axis.text.x = element_text(family="Arial", size = 12, color = "black"), #edit x-tick labels
        axis.text.y = element_text(family="Arial", size = 12, color = "black")) #I did it!
ggsave("loes_regression_MedNormFluo_011222.png")

### on HPC: For Loop - for each sample, subset it and write a sc_distributions_SampleName_allTimepoints.csv
#for(pop in unique(sc_distr_alltimepoints$sample)) {
#  print(pop)
#  sc_distr_alltimepoints %>%
#  filter(sample == pop) %>%
#  write_csv(paste0("sc_distributions_",pop,"_all_timepoints.csv"))
#}
#### Investigate the zig zaggy timepoints by graphing the ridgeplots for those populations to see what the distribution is like. (Zigzaggy lines of any one population can seen in the plot_list plots)
#my idea is maybe the distribution shape can tel whether it's CNV dynamics or contamination.
sc_gap1_4 = read.csv(file = "sc_distributions_gap1_4_all_timepoints.csv", stringsAsFactors = T) %>%
  mutate(generation = factor(generation, levels = unique(generation)))

sc_gap1_4 %>%
  ggplot(aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
  geom_density_ridges_gradient(scale = 2.0, rel_min_height = 0.01) +
  xlab("Normalized fluorescence (a.u.)") +
  ylab("Generation") +
  ggtitle("GAP1 WT 4") +
  theme_classic() +
  scale_x_continuous(limits=c(0.0,2.5), breaks = c(0, 1, 2, 2.5)) +
  scale_y_discrete(expand = expansion(add = c(0.2, 2.5))) + #expands the graph space or else the top is cut off
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
  theme(
    legend.position = 'none', #remove the legend
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
ggsave("gap1_4_ridgeplot_scale2.png")

sc_gap1_4 %>%
  group_by(generation) %>%
  mutate(norm_median = median(B2A_FSC)) %>%
  distinct() %>%
  ggplot(aes(generation, norm_median, group = 1)) +
  geom_line()+
  geom_point()+
  ggtitle("gap1_4")+
  ylab("normalized B2A/FSC then median")
  theme_classic()
ggsave("gap1_4_normalized_median_lineplot.png")

sc_gap1_4 %>%
  group_by(generation) %>%
  mutate(med_B2A = median(B2.A)) %>% View()
  ggplot(aes(generation, med_B2A, group = 1)) +
  geom_line() +
  geom_point()+
  ggtitle("gap1_4") +
  ylab("median B2-A fluorescence (a.u)") +
  theme_classic()
ggsave("gap1_4_raw-B2A_lineplot.png")

sc_gap1_4 %>%
  group_by(generation) %>%
  mutate(med_FSC = median(FSC.A)) %>%
ggplot(aes(generation, med_FSC, group = 1)) +
  geom_line() +
  geom_point()+
  ggtitle("gap1_4") +
  ylab("median FSC fluorescence (a.u)") +
  theme_classic()
ggsave("gap1_4_raw-FSC_lineplot.png")

clean_adj_norm_medians %>% filter(sample == "gap1_4") %>%
  ggplot(aes(generation, Med_B2A_FSC, group = 1))+
  geom_line()+
  geom_point()+
  theme_classic()
ggsave("gap1_4_median_then_norm_lineplot.png")

clean_adj_norm_medians %>%
  filter(sample == "gap1_4", generation > 120) %>%
  select(sample, Description, generation, `FSC-A`, GFP, Med_B2A_FSC, Count) %>% View()

#experimental populations that dip down, investigate.
dips = clean_adj_norm_medians %>%
  filter(Med_B2A_FSC < 1.5 & Type == "Experimental") %>%
  write_csv("Med_B2A_FSC_dips.csv")

#investigate sample gap1_all_3. Graph ridgeplots.
gap1_all_3 <- read.csv(file = "sc_distributions_gap1_all_3_all_timepoints.csv", stringsAsFactors = T) %>%
  mutate(generation = factor(generation, levels = unique(generation)))
gap1_all_3 %>%
ggplot(aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
  geom_density_ridges_gradient(scale = 2.0, rel_min_height = 0.01) +
  xlab("Normalized fluorescence (a.u.)") +
  ylab("Generation") +
  ggtitle("GAP1 LTR + ARS 3") +
  theme_classic() +
  #scale_x_continuous(limits=c(0.0,5), breaks = c(0, 1, 5, 0.5)) +
  scale_y_discrete(expand = expansion(add = c(0.2, 2.5))) + #expands the graph space or else the top is cut off
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
  theme(
    legend.position = 'none', #remove the legend
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
ggsave("gap1_all_3_ridgeplots.png")


####################################################
# STEP 9:  Quantify CNV dynamics (Lauer et al. 2018)
# Author: Julie
  # 1) First, calculate Tup, the generation at which CNVs are initially detected, (Lang et al. 2011 and Lauer et al. 2018)
  # To do that, calculate the false positive rate for CNV detection (threshold), which I will define as the median frequency of 1 copy control cells appearing in the two_copy_or_more gate and gate across generations 8-260 plus the interquartile range (IQR) if the distribution of one copy controls appearing in the CNV gate as NOT normal. If the distribution is normal, then use the mean plus one standard deviation like in Lauer et al. 2018. Like in Lauer et al. 2018, samples surpassing this threshold is considered to contain CNVs.

#CNV False Positive Rate is defined by the frequency of the 1 copy control strain appearing in the CNV gate which is called the Two_or_more copy gate.
CNV_false_pos_df = #freq %>%
  fw_freq_and_counts %>%
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
Tup_per_pop = #freq %>%
  fw_freq_and_counts %>%
  filter(Count>70000) %>%
  anti_join(fails) %>%
  filter(Type == "Experimental", Gate == "two_or_more_copy", Frequency >= thres_mean) %>%
  select(Type, Strain, Description, sample, generation, Gate, Frequency) %>% #View()
  group_by(sample) %>%
  slice(which.min(generation))
Tup_per_pop %>% write_csv("file = 01_02_04_v2_fw_Tup_per_pop.csv")

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
ggsave("01_02_04_v2_fw_Tup_boxplot.png")

# ANOVA to test for significance
# One Way ANOVA because there is only 1 independent variable, genotype.
# Anova Tut: https://www.scribbr.com/statistics/anova-in-r/
Tup_anova = aov(generation~Description, data = Tup_per_pop)
summary(Tup_anova)
#            Df  Sum Sq Mean Sq F value  Pr(>F)
#Description  3   4483  1494.3   6.772 0.00181 **
#Residuals   24   5296   220.7
# Conclusion: There IS a significant difference in the means. Genotype has has significant impact on Tup.

#Calculate Sup
#S1 Text. Calculation of CNV dynamics parameters.
# Graphic representation of linear fit (and corresponding R2 values)
# during initial population expansion of CNV alleles.
# Slope of the lineear fit corresponds to the dynamics parameter Sup shown in
# Table 1 and was calculated for the original evolution experiment and the
# barcode experiment. Data and code used to generate these figures can be
# accessed in OSF: https://osf.io/fxhze/. CNV, copy number variant.

#Use the "Generate S1 Text.Rmd" to calculate Sup
#I made my own wokring version - Generate S1 Text_JC.Rmd inside the cnvflow folder
# from current freq table, I need to add the following 4 columns:
# propCNV, propNoCNV, propCNV_divided_by_propNoCNV, natural log of propCNV_divided_by_propNoCNV
# not log base10, remember that natural log is ln base e. ln x= log base e of x.
# log makes it linear -- and we are doing a linear fit
# why use natural log? I don't think there's a biological reason. I think the reason is mathematical,
# in that the derivative or slope of an natural log line y=ln(x) is simply 1/x. In other words dy/dx ln(x) = 1/x.

# Step 1 - make a table like Steff's FlowAnalysisSumm_FINAL.csv. For me that's the

equation = function(x) {
  lm_coef <- list(a = round(coef(x)[1], digits = 2),
                  #b = round(coef(x)[2], digits = 2),
                  b = unname(round(coef(fit)[2], digits = 2)),# get rid of it printing c()
                  r2 = round(summary(x)$r.squared, digits = 2));
  lm_eq <- substitute(slope == b~~~~italic(R)^2~"="~r2,lm_coef)
  as.character(as.expression(lm_eq));
}
#format(unname(coef(m))[1], digits = 2)
dynamics = fw_freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(fails)  %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  group_by(sample, generation) %>%
  mutate(prop_CNV = sum(Frequency),
         prop_NoCNV = 100-prop_CNV,
         CNV_NoCNV = prop_CNV/prop_NoCNV,
         logECNV_NoCNV = log(CNV_NoCNV)) #log() function is natural logarithm in R (even though  log() commonly thought as base10 )

#write a function to calculate Sup, Explained Variance, make graphs, ggsave graphs
#then, use map() to apply function to all populations - I have 28
#the tricky thing is the generations bounds can be different for each population
pop_list = unique(dynamics$sample)
wt_pops = pop_list[c(25,21,26,22,27)] #subset the list
pop_data <- subset(dynamics, sample %in% c(pop_list[[27]]) & generation >=29 & generation <=124) #why did steff choose gen 41 - gen124?
fit <- lm(logECNV_NoCNV ~ generation, pop_data) #linear model, lm(y~x,data)
fit
#summary(fit) #to see the full model
ggplot(subset(dynamics, sample %in% c(pop_list[[27]])), aes(x=generation,y=(as.numeric(logECNV_NoCNV)), colour=sample)) +
  geom_point() +
#  geom_smooth(data=subset(pop_data, generation >=41 & generation <=124), method=lm, show.legend=FALSE) +
  geom_smooth(data=pop_data, method=lm, show.legend=FALSE) +
 # scale_y_continuous(expand = c(0, 0), 'ln(Prop. CNV/Prop. non-CNV)', limits=c(-5,3)) +
  scale_y_continuous(expand = c(0, 0), 'ln(Prop. CNV/Prop. non-CNV)', limits = c(min(pop_data$logECNV_NoCNV)-1, max(pop_data$logECNV_NoCNV)+1)) +
  annotate("text", x = 130, y = 0.5, label = equation(fit), parse = TRUE) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5), "Generations", limits=c(0,260)) +
  theme_classic() +
  scale_color_manual(values = c('black')) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  theme(legend.position = c(.15,.95), plot.title = element_text(size=14, hjust = 0.5), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))

summary(fit)$coef[[4]]

confint(fit, 'Generation', level = 0.95)





####### MY PALLETEE
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

