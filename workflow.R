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
library(docstring)

# Set working directory and get list of subdirectories containing FCS files
#setwd('/Volumes/GoogleDrive/My Drive/Gresham Lab_Papers/2021/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files') #David's working directory
#setwd('G:/.shortcut-targets-by-id/1Bioj1YP_I7P8tqgmg4Zbt4EAfhb7J-0w/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files') #Titir's working directory
#setwd("/Volumes/GoogleDrive/My Drive/greshamlab/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files") #Julie's WD
setwd("/Volumes/GoogleDrive/My Drive/greshamlab/projects/EE_GAP1_ArchMuts_Summer2021/data/Summer_LTEE_2021_FCS_files")  #Julie's WD

#In addition to having  directories containing data FSC files, make a gating directory, which is a directory that contains ALL the FSC files you want to overlay for drawing gates.
folders = list.dirs()[c(2,5:28)] #select the FSC file folders in your directory

# Choose a name to be used for all output files including the gating template and associated flow data and graphs.
version_name = "01_02_04_v4_wt_ltr" #timepoints folders 1,2,4 using WT and LTR and control samples only to draw gates (folder 01_02_04_wt_ltr_ctrls_only_FSCfiles) using 0,1,2 copy as guides.
#version_name = "01_02_04_v3_wt_ltr" #timepoints folders 1,2,4 [using WT and LTR and control samples only](<- i don't think this is true. I think all samples were used to draw gates) using 0,1 copy as guides.
#version_name = "newGates_01_02_04_ars_all" #timepoint tolders 1,2,4, using only ARS and LTR+ARS samples (folder 01_02_04_ars_all_only_FSCfiles) only to draw gates
#version_name = "01_02_04_v2" #timepponts folders 1,2,4 all samples (experimental and controls) to draw gates
# other versions: 01_02_04_v2

#STEP 1: Generate experiment details file from folder and FCS file names
# Experiment details file is a .csv file that contains the list of .fcs files in the directory and the associated metadata for each sample
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

#needs to be run once
map(folders, make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021.csv")


#STEP 2: Read in all files in a directory and rename the channels.
#A directory contains an FCS file for each population.
#Results in 1 timepoint gating set containing all .fcs files, associated experiment details, and marker details
#Author: Julie

# here we will load in 1 directory that contains 3 timepoints worth of data to load in.
# cyto_setup() does not permit loading in more than 1 directory, so I had to create a directory with the data files of interest.
# these data will guide us on drawing gates.
exp_details_path = list.files(path = paste0(folders[1]), pattern = "_experiment_details.csv", full.names = T)

timepoint_gating_set <- cyto_setup(path = paste0(folders[1]), restrict=TRUE, select="fcs", details=F) #edit Markers on Viewer pane, Save & Close

#use flowWorkspace::pData to annotate the experiment details file associated with the gating set
experiment_details <- read_csv(exp_details_path) #import experiment-details.csv

for(i in 1:length(names(experiment_details))){
  flowWorkspace::pData(timepoint_gating_set)[names(experiment_details[i])]<-experiment_details[i]
  }

## Rename the experiment-markers.csv file. Need to do once.
#file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv")

#STEP 3:  Perform gating on gating set
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.
#Author: Titir & Julie

## 3.1 transform the data
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

## quickly check the transformation by plotting the data
#my_samples = which(experiment_details$Description %in% c("0 copy control","1 copy control","2 copy control","GAP1 LTR KO", "GAP1 WT architecture"))
cyto_plot_explore(transformed_timepoint_gating_set,
                  channels_x = "FSC-A",
                  channels_y = "B2-A",
                  axes_limits = "data")

## 3.2 Gating using the entire timepoint dataset or apply an existing gating template

# note:if you already have a gating template and don't need to draw gates, then skip cyto_draw, use cyto_gatingTemplate_apply to apply the gating template.csv to your gating set
cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= "cytek_gating_01_02_04_v2.csv")

#First we gate for the cells
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = paste0("cytek_gating_",version_name,".csv")
)


#Then we define the singlets based on forward scatter height and width
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = paste0("cytek_gating_",version_name,".csv")
)

#Gating for CNVs using the 0,1 and 2 copy controls:
indexes_ctr0 <- which(experiment_details$Description %in% c("0 copy control"))
DGY1 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[as.numeric(indexes_ctr0)] #DGY1 c(30,61,92)

ind_ctr1 <-as.numeric(which(experiment_details$Description %in% c("1 copy control")))
DGY500 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[ind_ctr1] #DGY500

ind_wt_ltr <-as.numeric(which(experiment_details$Description %in% c("GAP1 LTR KO", "GAP1 WT architecture")))
exp_1_copy <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[ind_wt_ltr] # experimental strains which I assume are 1 copy.

indexes_ctr2 <- as.numeric(which(experiment_details$Description %in% c("2 copy control")))
DGY1315 <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[indexes_ctr2] #DGY1315


cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Single_cells", #first color
               alias = c("zero_copy", "one_copy", "two_or_more_copy"), #defines gate names
               channels = c("FSC-A","B2-A"),
               axes_limits = "data",
               #select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
               gatingTemplate = paste0("cytek_gating_",version_name,".csv"),
#               overlay = c(DGY1, DGY500, DGY1315),
               overlay = c(DGY1, DGY500, exp_1_copy, DGY1315),
               point_col = c("gray", "green", "black", "purple","blue") #parent color then overlay colors
)

#STEP 4:  Generate single cell data tables and normalized fluorescence
#Author: Julie

timepoint_raw_list <- cyto_extract(transformed_timepoint_gating_set, parent = "Single_cells", raw = TRUE, channels = c("FSC-A", "B2-A")) #raw flow data of each single cell as a list of matrices

map_df(timepoint_raw_list, ~as.data.frame(.x), .id="name") %>% #convert to df, put list name in new column
  mutate(name = as.factor(name)) %>% #convert `name` to factor
  left_join(experiment_details %>% #join by name column to add metadata
  mutate(generation = as.factor(unique(experiment_details$generation)))) %>%
  mutate(B2A_FSC = `B2-A`/`FSC-A`) %>% #compute normalized fluor
  write_csv(paste0(version_name,"_SingleCellDistributions_",prefix,".csv"))

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

  #6. Get cell counts and frequencies inside each gate
  #Julie

  #get cell count from each gate
  gs_pop_get_stats(transformed_timepoint_gating_set, c("Single_cells", "zero_copy", "one_copy", "two_or_more_copy")) %>%
    rename(Gate = pop, name = sample, Count = count) %>%
    left_join(experiment_details) %>%
    write_csv(paste0(version_name,"_counts_",prefix,".csv"))

  #get frequency of cells inside each gate
  gs_pop_get_stats(transformed_timepoint_gating_set, c("Single_cells","zero_copy", "one_copy", "two_or_more_copy"), type = "percent") %>%
    rename(Gate = pop, name = sample, Frequency = percent) %>%
    left_join(experiment_details) %>%
    write_csv(paste0(version_name,"_freq_",prefix,".csv"))

  #get single cell fluorescence normalized over cell size
#  timepoint_raw_list <- cyto_extract(transformed_timepoint_gating_set, parent = "Single_cells", raw = TRUE, channels = c("FSC-A", "B2-A")) #raw flow data of each single cell as a list of matrices

#  map_df(timepoint_raw_list, ~as.data.frame(.x), .id="name") %>% #convert to df, put list name in new column
#    mutate(name = as.factor(name)) %>% #convert `name` to factor
#    left_join(experiment_details %>% #join by name column to add metadata
#                mutate(generation = as.factor(unique(experiment_details$generation)))) %>%
#    mutate(B2A_FSC = `B2-A`/`FSC-A`) %>% #compute normalized fluor
#    write_csv(paste0(version_name,"_SingleCellDistributions_",prefix,".csv"))


}

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Julie

try(map(folders[2:length(folders)],analyze_all_exp, my_markers, gating_template = paste0("cytek_gating_",version_name,".csv")))

#STEP 7: Pull in all counts or freq or single cell distribution files from directory and combine into a single dataframe
#Author: Julie
#01_02_04_v2_fw_counts_all_timepoints
#01_02_04_v2_fw_freq_all_timepoints
list.files(path = ".", pattern = paste0(version_name,"_counts_([0-9])+_EE_GAP1_ArchMuts_2021")) %>%
  read_csv() %>%
  mutate(gating_template = paste0("cytek_gating_",version_name,".csv")) %>%
  write_csv(file = paste0(version_name,"_counts_all_timepoints.csv")) #write a function to establish/customize the version file name then can paste0(version_filename)

list.files(path = ".", pattern = paste0(version_name,"_freq_([0-9])+_EE_GAP1_ArchMuts_2021")) %>%
  read_csv() %>%
  mutate(gating_template = paste0("cytek_gating_",version_name,".csv")) %>%
  write_csv(file = paste0(version_name,"_freq_all_timepoints.csv"))

## Do on hpc because large files, do once. Don't even try to run this command on your laptop. 12GB file.
# list.files(path = ".", pattern = paste0(version_name,"_SingleCellDistributions")) %>%
#   read_csv() %>%
#   mutate(gating_template = paste0("cytek_gating_",version_name,".csv")) %>%
#   write_csv(file = paste0(version_name,"_SingleCellDistributions_all_timepoints.csv"))

#STEP 8: Plot cells in gates ridgeplots, time series, & assess gates
#Determine whether =>83% of controls are in the correct gate
#Make plots
#Author: Grace & Julie

# read in frequency csv, cell numbers csvs, single cell distributions for all timepoints
#freq = read_csv(paste0(version_name,"_freq_all_timepoints.csv"))
#freq = read_csv("01_02_04_v2_fw_freq_all_timepoints.csv")
#freq = read_csv("newGates_01_02_04_ars_all_freq_all_timepoints.csv")


#count= read_csv(paste0(version_name,"_counts_all_timepoints.csv"))
#count= read_csv("01_02_04_v2_fw_counts_all_timepoints.csv")
#count=read_csv("newGates_01_02_04_ars_all_counts_all_timepoints.csv")

freq_and_counts = read_csv("freq_and_counts_Merged_080622_all_timepoints.csv")

sc_distr_alltimepoints <- read.csv(paste0(version_name,"_SingleCellDistributions_all_timepoints.csv", stringsAsFactors = T)) %>% mutate(generation = factor(generation, levels = unique(generation)))

freq_and_counts =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate) %>%
  relocate(9, .after = Frequency)

#Table of low cell observations, convenient to have to anti_join() in further steps
lowcell = freq_and_counts %>%
  filter(Count <7000) %>%
  mutate(generation = factor(generation, levels = unique(generation))) %>% #View()
  select(-Count)

## check controls are in their proper gates
  fails = freq_and_counts %>%
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
  dplyr::filter(flag == "fail") %>%
  arrange(Description)
  View(fails)
  #fails %>% write_csv("01_02_04_v2_83_fail.csv")
  #fails %>% write_csv("01_02_04_v2_fail_calc_thres_stringent_.csv")
  #fails %>% write_csv("01_02_04_v2_79_10_fail_.csv")
  #fails %>% write_csv("01_02_04_v2_fw_79_11_fail.csv")

# plot proportion of control cells in control gates over time
freq_and_counts %>%
  filter(Count>70000,
          str_detect(Description, "control"),
         generation <250) %>%
  select(Type, Strain, Description, generation, Gate, Frequency, Count) %>%
  dplyr::filter(!(Description == "1 copy control" & generation == 182 |
                  Description == "2 copy control" & generation == 79 |
                  Description == "2 copy control" & generation == 95 |
                  Description == "2 copy control" & generation == 108 |
                  Description == "2 copy control" & generation == 116)) %>% #exclude these controls timepoints that look weird on ridgeplots
  #anti_join(fails) %>% #exclude the contaminated controls timepoints (the failed timepoints)
  ggplot(aes(generation, Frequency, color = Gate)) +
  geom_line() +
  facet_wrap(~Description) +
  ylab("% of cells in gate") +
  theme_minimal() +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme(text = element_text(size=12))

ggsave(paste0("propCNV_",version_name,"_controls_8x12.pdf"), bg = "#FFFFFF", height = 8, width = 12)
ggsave(paste0("propCNV_",version_name,"_controls_10x14.pdf"), bg = "#FFFFFF", height = 10, width = 14)

# plot proportion of population in each gate over time for each of 28 experimental populations
prop_plot_list = list()
i=1
for(exp in unique(freq_and_counts$Description)) {
  prop_plot_list[[i]] = freq_and_counts %>%
    filter(Count>70000) %>%
    #filter(generation != 79, generation != 116,generation != 182,generation != 252) %>%
    filter(Description==exp) %>%
    ggplot(aes(generation, Frequency, color = Gate)) +
    geom_line(size =1.5) +
    facet_wrap(~sample) +
    ylab("% of cells in gate") +
    theme_minimal()+
    scale_x_continuous(breaks=seq(0,250,50))+
    theme(text = element_text(size=10))
  i = i+1
}
names(prop_plot_list) = unique(freq_and_counts$Description)
prop_plot_list$`GAP1 WT architecture` # change index to view replicates for different genetic backgrounds
prop_plot_list$`GAP1 ARS KO`
prop_plot_list$`GAP1 LTR KO`
prop_plot_list$`GAP1 LTR + ARS KO`


ggsave(paste0("propCNV_",version_name,"_ALL.pdf"), bg = "#FFFFFF", height = 5, width = 12)
### Plot proportion of the population with a CNV over time

my_facet_names <- as_labeller(c("GAP1 WT architecture" = "Wildtype architecture",
                          "GAP1 LTR KO" = "LTR KO",
                        "GAP1 ARS KO" = "ARS KO",
                        "GAP1 LTR + ARS KO" = "LTR and ARS KO"))
#colors
# wtGrays = c("gray","#666666","#CCCCCC","gray","#999999")  #OLD
wtGrays = c("#354f52","#666666","#6b705c","#414833","#999999")
#"#354f52", "#414833","#6b705c"
# allGolds = c("#DEBD52","#DBB741","#D7B02F","#dbb844","#D9BB59","#fdc409","#9c7e1e","#D9BB59") #OLD
allGolds = c("#ffba08", "#faa307", "#dda15e", "#7f5539", "#9c6644", "#fdc409", "#9c7e1e","#D9BB59")
#"#dda15e" #nude
#"#e85d04" orange-red
#ee9b00 #gold
#ca6702 #pumpkin
#bb3e03 #warmer pumpkin
#ae2012 #dark red
#keep fdc409 its the super bright yellow one
#keep 9c7e1e its the brown one
#add some browns to the yellows  #"#b08968", "#7f5539", "#9c6644"
arsSalmons = c("#e26d5c","#e28f5c","#e25c6d","#da4631", "#f85c46", "#bb3521","#d9402a" )
#ltrBlues = c("#6699cc", '#66b3cc',"#6BAED6" ,"#4292C6", "#2171B5","#3799fb","#3972ab","#4799eb") #old

ltrBlues = c( "#6699cc", "#005f73", "#0a9396", "#4292C6", "#2171B5", "#3799fb", '#66b3cc', "#3a0ca3")

propCNV = freq_and_counts %>%
  filter(Count>70000,
  generation <= 203) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(fails)  %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  dplyr::filter(!(Description == "1 copy control" & generation == 182 |
                    Description == "2 copy control" & generation == 79 |
                    Description == "2 copy control" & generation == 95 |
                    Description == "2 copy control" & generation == 108 |
                    Description == "2 copy control" & generation == 116)) %>% #exclude these controls timepoints that look weird on ridgeplots
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_line(size = 2.5) +
  #geom_point()+
  facet_wrap(~factor(Description,
              levels = c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")), labeller = my_facet_names, scales='free') +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications") +
  scale_color_manual(values = c(wtGrays, allGolds,arsSalmons, ltrBlues)) +
  theme_classic() +
  #scale_x_continuous(breaks=seq(0,250,50)) +
  scale_x_continuous(breaks=seq(0,203,50)) +
  scale_y_continuous(limits=c(0,100)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size=25),
        legend.position = "none",
        axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=25)
        )
propCNV

ggsave(paste0("propCNV_",version_name,"_080722_8x12.pdf"), bg = "#FFFFFF", height = 8, width = 12)
ggsave(paste0("propCNV_",version_name,"_080722_10x14.pdf"), bg = "#FFFFFF", height = 10, width = 14)
ggsave("propCNV_101322_8x12.pdf", bg = "#FFFFFF", height = 8, width = 12)
ggsave("propCNV_101322_8x12.png", bg = "#FFFFFF", height = 8, width = 12)
ggsave("propCNV_101322_10x14.pdf", bg = "#FFFFFF", height = 10, width = 14)


# proportion of population in CNV over time for each population
# facet by sample instead of Description
propCNV_by_Pop = freq_and_counts %>%
  filter(Count>70000,
         # generation <= 203) %>%
         generation <= 203) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(fails)  %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  dplyr::filter(!(Description == "1 copy control" & generation == 182 |
                    Description == "2 copy control" & generation == 79 |
                    Description == "2 copy control" & generation == 95 |
                    Description == "2 copy control" & generation == 108 |
                    Description == "2 copy control" & generation == 116)) %>% #exclude these controls timepoints that look weird on ridgeplots
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_line(size = 2.5) +
  #geom_point()+
  facet_wrap(~sample, scales='free') +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications") +
  scale_color_manual(values = c(
    wtGrays,
    allGolds,
   # "#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59","#D7B02F","#CAA426","#D9BB59", #ALL ko ,8,gold
arsSalmons,
    ltrBlues
  )) +
  theme_minimal() +
  scale_x_continuous(breaks=seq(0,203,50)) +
  scale_y_continuous(limits=c(0,100)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size=12),
        legend.position = "none",
        axis.text.x = element_text(size = 10, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 12, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=12)
  )
propCNV_by_Pop

ggsave(paste0("propCNV_by_pop_",version_name,"_080722.pdf"), bg = "#FFFFFF", height = 15, width = 20)
ggsave("propCNV_by_pop_101322.pdf", bg = "#FFFFFF", height = 15, width = 20)
ggsave("propCNV_by_pop_101322.png", bg = "#FFFFFF", height = 15, width = 20)


############
# PropCNV plots with abberant timepoints and/or populations removed

early_df = freq_and_counts %>%
  filter(generation < 30,
         Type %in% c("Experimental", "1_copy_ctrl"),
         Description %in% c("1 copy control", "GAP1 WT architecture","GAP1 LTR KO"),
         Gate == "two_or_more_copy") %>%
  arrange(generation, sample)%>%
  select(-name, -`Outflow well`, -Media)

weird_early = early_df %>%
  filter(Frequency > 15)

freq_and_counts %>%
  filter(sample == "gap1_all_6", Gate == "two_or_more_copy") %>%
  arrange(generation) %>%
  View()

#chose these timepoints by eye
weird_tp = freq_and_counts %>%
  filter(sample == "gap1_4" & Gate == "two_or_more_copy" & generation == 66 |
        sample == "gap1_all_3" & Gate == "two_or_more_copy" & generation == 166|
        sample == "gap1_all_5" & Gate == "two_or_more_copy" & generation == 116|
        sample == "gap1_all_6" & Gate == "two_or_more_copy" & generation == 124|
        sample == "gap1_ltr_2"
    )

freq_and_counts %>%
  filter(Count>70000,
         generation <= 203) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(fails)  %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  dplyr::filter(!(Description == "1 copy control" & generation == 182 |
                    Description == "2 copy control" & generation == 79 |
                    Description == "2 copy control" & generation == 95 |
                    Description == "2 copy control" & generation == 108 |
                    Description == "2 copy control" & generation == 116)) %>% #exclude these controls timepoints that look weird on ridgeplots
  anti_join(weird_early) %>%
  anti_join(weird_tp) %>%
  ggplot(aes(generation, Frequency, color = sample)) +
  geom_line(size = 2.5) + #alpha=0.7
  #geom_point()+
  facet_wrap(~factor(Description,
                     levels = c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")), labeller = my_facet_names, scales='free') +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplifications") +
  scale_color_manual(values = c(wtGrays, allGolds,arsSalmons, ltrBlues)) +
  theme_classic() +
  #scale_x_continuous(breaks=seq(0,250,50)) +
  scale_x_continuous(breaks=seq(0,203,50)) +
  scale_y_continuous(limits=c(0,100)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size=25),
        legend.position = "none",
        axis.text.x = element_text(size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(size = 30, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=25)
  )




##### STEP 9 ##### Plot Ridgeplots (density histograms):
# instead of looking at the MEDIAN GFP values per population per generation,
# we want to see the ENTIRE DISTRIBUTION of GFP of single cells per population per generation.

#Plot ridgeplots of each experimental and control population (32 in all) to catch contaminated timepoints/outlier values.
# Wrote a Function that plots ridgeplots given the single_cell_distribution_{POP-NAME}_all_timepoints.csv
# Then I can write use map() to apply this FUNCtion to all population.csv since I have 32 pops

pop_files = list.files(pattern = "sc_distributions_") %>% sort()
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

map(pop_files, make_ridgeplots) #map() applies this ridgeplot function to all 32 population.csv files

### on HPC: For Loop - for each sample, subset it and write a sc_distributions_SampleName_allTimepoints.csv
#for(pop in unique(sc_distr_alltimepoints$sample)) {
#  print(pop)
#  sc_distr_alltimepoints %>%
#  filter(sample == pop) %>%
#  write_csv(paste0("sc_distributions_",pop,"_all_timepoints.csv"))
#}



#Make Faceted Ridgeplots of Control strains over time like Fig2A in Lauer et al. 2018.

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


### STEP 10
###### Plot Median normalized GFP fluorescence plot with controls separate
## 5 plots

# on hpc, do once
sc_distr_alltimepoints %>%
  group_by(sample, generation) %>%
  mutate(Med_B2A_FSC = median(B2A_FSC)) %>%
  distinct(Med_B2A_FSC, .keep_all = T) %>%
  select(-FSC.A, -B2.A, -B2A_FSC) %>%
  write_csv("medians_normalized_fluor_alltimepoints.csv")

cell_numbers = count %>%
  filter(Gate == "Single_cells")

norm_medians = read_csv("medians_normalized_fluor_alltimepoints.csv") %>%
  left_join(cell_numbers) %>%
  filter(Count > 70000) %>% #filter out low cells
  filter(!(generation == 231 & Type == "0_copy_ctrl")) %>%   #filter out bad controls
  filter(!(generation == 182 & Type == "1_copy_ctrl")) %>%
  filter(!(generation == 203 & Type == "1_copy_ctrl")) %>%
  filter(!(generation == 252 & Type == "1_copy_ctrl")) %>%
  filter(!(generation == 260 & Type == "1_copy_ctrl")) %>%
  filter(!(generation == 79 & Type == "2_copy_ctrl")) %>%
  filter(!(generation == 95 & Type == "2_copy_ctrl")) %>%
  filter(!(generation == 108 & Type == "2_copy_ctrl")) %>%
  filter(!(generation == 116 & Type == "2_copy_ctrl"))


MedianFluor_plot = norm_medians %>%
  filter(Type == "Experimental") %>%
  filter(!(Med_B2A_FSC<1.5 & Type == "Experimental")) %>%  #filter out outliers (likely resulting from contamination) as defined by Fluor <1.5
  ggplot(aes(generation, Med_B2A_FSC, color= sample)) +
  geom_line(size = 2.0) +
  scale_color_manual(values = c(wtGrays, allGolds, arsSalmons, ltrBlues)) +
  facet_wrap(~factor(Description,
                     levels = c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")), labeller = my_facet_names, scales='free') +
  xlab("Generation") +
  ylab("Median normalized fluorescence (a.u.)") +
  scale_x_continuous(breaks=seq(0,200,50)) +
  xlim(0,225)+
  ylim(c(1.5,2.5))+
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=36),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=36),
        #axis.text.x = element_text(family="Arial", size = 24, color = "black"), #edit x-tick labels
        axis.text.x = element_text(size = 36, color = "black"), #edit x-tick labels
        #axis.text.y = element_text(family="Arial", size = 24, color = "black")
        axis.text.y = element_text(size = 36, color = "black")
  )
MedianFluor_plot
ggsave("medNormFluorPlot.pdf")
ggsave("medNormFluorPlot.png")


### graph the controls separately
MedianFluor_ctrl = norm_medians %>%
  filter(Type %in% c("2_copy_ctrl","1_copy_ctrl")) %>%
  #filter(!(Med_B2A_FSC<1.5 & Type == "Experimental")) %>%  #filter out outliers (likely resulting from contamination) as defined by Fluor <1.5
  ggplot(aes(generation, Med_B2A_FSC, color= sample)) +
  geom_line(size = 2.0) +
  scale_color_manual(values = c("black", "black"))+
  xlab("Generation") +
  ylab("Median normalized fluorescence (a.u.)") +
  scale_x_continuous(breaks=seq(0,200,50)) +
  xlim(0,225)+
  ylim(c(1.5,2.5))+
  theme_classic() +
  ggtitle("1 and 2 copy controls") +
  theme(legend.position = "none",
        #legend.title = element_blank(),
        text = element_text(size=36),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=36),
        #axis.text.x = element_text(family="Arial", size = 24, color = "black"), #edit x-tick labels
        axis.text.x = element_text(size = 36, color = "black"), #edit x-tick labels
        #axis.text.y = element_text(family="Arial", size = 24, color = "black")
        axis.text.y = element_text(size = 36, color = "black")
  )
MedianFluor_ctrl
ggsave("controls_normMedFluorPlot.png")


# Graph single plots of median normalized fluorescence for every population (sample)
# Later, match it up with the single cell distribution ridgeplots to see if ridgeplots are consistent with the median.
# (They should be consistent)
fluor_single_plots = list()
i=1
for(exp in unique(clean_adj_norm_medians$Description)) {
  fluor_single_plots[[i]] = clean_adj_norm_medians %>%
    #filter(Description==exp) %>%
    ggplot(aes(generation, Med_B2A_FSC, color= sample)) +
    geom_line(aes(linetype = Type), size = 1.5) +
    facet_wrap(~sample) +
    scale_linetype_manual(values = c("dashed", "dashed", "dashed", "solid")) +
    xlab("Generation") +
    ylab("Median normalized fluorescence (a.u.)") +
    scale_x_continuous(breaks=seq(0,250,50))+
    theme_bw() +
    theme(#text = element_text(size=20),
          axis.text.x = element_text(size=10),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          strip.background = element_blank(), #removed box around facet title
    )
  i = i+1
}
names(fluor_single_plots) = unique(clean_adj_norm_medians$Description)
fluor_single_plots$`GAP1 WT architecture` # change index to view replicates for different genetic backgrounds
fluor_single_plots$`GAP1 ARS KO`
fluor_single_plots$`GAP1 LTR KO`
fluor_single_plots$`GAP1 LTR + ARS KO`


####################################################
# STEP 11:  Quantify CNV dynamics (Lauer et al. 2018)
# see script quant_cnv_dynamics.R
# Author: Julie

###########  MY PALLETTE #####

# for lineplots population-specific

#Gold Metallic (6)
#DEBD52
#DBB741
#D7B02F
#CAA426
#D9BB59
#B89523

#colors for lineplots 08-07-22
# wtGrays = c("gray","#666666","#CCCCCC","gray","#999999")
# allGolds = c("#DEBD52","#DBB741","#D7B02F","#dbb844","#D9BB59","#fdc409","#9c7e1e","#D9BB59")
# arsSalmons = c("#e26d5c","#e28f5c","#e25c6d","#da4631", "#f85c46", "#bb3521","#d9402a" )
# ltrBlues = c("#6699cc", '#66b3cc',"#6BAED6" ,"#4292C6", "#2171B5","#3799fb","#3972ab","#4799eb")

#NEW PALLETE 4/20/22 use this for boxplots
#WT = gray  "gray","gray","gray","gray","gray",
#LTR and ARS = GOLD metalic "#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59","#D7B02F","#CAA426","#D9BB59", #LTR,8,gold
#ARS = SALMON  "#e26d5c"
#LTR = BABY BLUE "#6699cc"
