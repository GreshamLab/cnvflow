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

folders = list.dirs()[5:28] #select the FSC file folders in your directory

# Choose a name to be used for all output files including the gating template and associated flow data and graphs.
version_name = "newGates_080322"

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
#A directory contains an FCS file for each populations for a single timepoint.
#Results in 1 timepoint gating set containing all .fcs files, associated experiment details, and marker details
#Author: Julie

# here we will choose 3 timepoints worth of data to load in
# these data will guide us drawing gates.
exp_details_path = list.files(path = paste0(folders), pattern = "_experiment_details.csv", full.names = T)

timepoint_gating_set <- cyto_setup(path = paste0(folders[1]), restrict=TRUE, select="fcs", details=F) #edit Markers on Viewer pane, Save & Close

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
               gatingTemplate = "cytek_gating_",version_name,".csv"
)


#Then we define the singlets based on forward scatter height and width
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "cytek_gating_",version_name,".csv"
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
               gatingTemplate = "cytek_gating_",version_name,".csv",
               overlay = c(DGY1, DGY500, DGY1315),
               point_col = c("gray", "green", "red", "blue")
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
}

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Julie

try(map(folders[1:length(folders)],analyze_all_exp, my_markers, gating_template = "cytek_gating_01_02_04_v2.csv"))

#STEP 7: Pull in all counts or freq or single cell distribution files from directory and combine into a single dataframe
#Author: Julie

list.files(path = ".", pattern = paste0(version_name,"_counts_([0-9])+_EE_GAP1_ArchMuts_2021")) %>%
  read_csv() %>%
  write_csv(file = paste0(version_name,"_counts_all_timepoints.csv")) #write a function to establish/customize the version file name then can paste0(version_filename)

list.files(path = ".", pattern = paste0(version_name,"_freq_([0-9])+_EE_GAP1_ArchMuts_2021")) %>%
  read_csv() %>%
  write_csv(file = paste0(version_name,"_freq_all_timepoints.csv"))

# Do on hpc because large files, do once
# list.files(path = ".", pattern = "01_02_04_v2_SingleCellDistributions") %>%
#   read_csv() %>%
#   write_csv(file = "01_02_04_v2_SingleCellDistributions_all_timepoints.csv")

#STEP 8: Plot ridgeplots, time series, & assess gates
#Determine whether =>83% of controls are in the correct gate
#Make plots
#Author: Grace & Julie

# read in frequency csv, cell numbers csvs, single cell distributions for all timepoints
freq = read_csv("01_02_04_v2_fw_freq_all_timepoints.csv") #%>% rename(Frequency = frequency)
count= read_csv("01_02_04_v2_fw_counts_all_timepoints.csv")
sc_distr_alltimepoints <- read.csv("01_02_04_v2_SingleCellDistributions_all_timepoints.csv", stringsAsFactors = T) %>% mutate(generation = factor(generation, levels = unique(generation)))
freq_and_counts =
  count %>% filter(Gate == "Single_cells") %>%
  rename(Parent = Gate) %>%
  left_join(freq) %>%
  filter(!(Gate == "Single_cells")) %>%
  mutate(Frequency = Frequency*100) %>%
  relocate(2:3, .after = Gate)

#Table of low cell observations, convenient to have to anti_join() in further steps
lowcell = freq_and_counts %>%
  filter(Count <7000) %>%
  mutate(generation = factor(generation, levels = unique(generation))) %>% #View()
  select(-Count)

## check controls are in their proper gates
  fails = freq_and_counts %>%
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

# plot proportion of control cells in control gates over time
freq_and_counts %>%
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


# Plot proportion of the population with a CNV over time

my_facet_names <- as_labeller(c("GAP1 WT architecture" = "Wildtype architecture",
                          "GAP1 LTR KO" = "LTR KO",
                        "GAP1 ARS KO" = "ARS KO",
                        "GAP1 LTR + ARS KO" = "LTR and ARS KO"))

freq_and_counts %>%
  filter(Count>70000,
         generation <= 203) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(fails)  %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  group_by(sample, generation) %>%
  mutate(prop_CNV = sum(Frequency)
         ) %>% #View()
  select(sample, generation, Description, prop_CNV) %>%
  distinct() %>%
  ggplot(aes(generation, prop_CNV, color = sample)) +
  geom_line(size = 2.5) +
  #geom_point()+
  facet_wrap(~factor(Description,
              levels = c("GAP1 WT architecture","GAP1 LTR KO", "GAP1 ARS KO","GAP1 LTR + ARS KO")), labeller = my_facet_names, scales='free') +
  xlab("Generation") +
  ylab("Proportion of cells with GAP1 amplications") +
#  scale_color_manual(values = c("#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59", #WT,5,gold
#"#637EE7","#6F88E9","#7B92EA","#4463E2","#3053DF","#2246D7","#1E3FC3","#5766E6", #ALL,8,bluepurple
#"#DE54B9","#E160BE","#E36CC3","#E578C8","#E885CD","#DB41B2","#D72FAA", #ARS,7,pink
#"#54DE79","#41DB6A","#2FD75C","#26CA52","#23B84B","#60E182","#6CE38C","#78E595" #LTR,8,green
#   )) +
  scale_color_manual(values = c(
  "gray","gray","gray","gray","gray", #wildtype, 5, gray
  "#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59","#D7B02F","#CAA426","#D9BB59", #LTR,8,gold
  "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", #ARS, 7, softer salmon repeats
  "#6699cc","#6699cc","#6699cc","#6699cc","#6699cc","#6699cc","#6699cc","#6699cc" #LTR,8,
)) +
  theme_classic() +
  scale_x_continuous(breaks=seq(0,250,50)) +
  scale_y_continuous(limits=c(0,100)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        text = element_text(size=25),
        legend.position = "none",
        axis.text.x = element_text(family="Arial", size = 30, color = "black"), #edit x-tick labels
        axis.text.y = element_text(family="Arial", size = 30, color = "black"),
        strip.background = element_blank(), #removed box around facet title
        strip.text = element_text(size=25)
        )

###################################################
# Plot Ridgeplots (density histograms):
# instead of looking at the MEDIAN GFP values per population per generation,
# we want to see the ENTIRE DISTRIBUTION of GFP of single cells per population per generation.

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

#Plot ridgeplots of each experimental and control population (32 in all) to catch contaminated timepoints/outlier values.
#Wrote a Function that plots ridgeplots.
# Then I can write use map() to apply this FUNCtion to all population.csv since I have 32 pops

pop_files = list.files(pattern = "sc_distributions_")[-1:-3] %>% sort()
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

cell_numbers = count %>%
  filter(Gate == "Single_cells")

norm_medians = read_csv("medians_normalized_fluor_alltimepoints.csv") %>%
  left_join(cell_numbers)
#%>%
  #select(-Marker) %>%
  #rename(Gate = Population)

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
#anti_join() gap1_4 two copy gate samples at g79, 124, 231, 252 -- why? not justified. That's what the data show. We don't know whether it's real or from contamination.

#Graph experimental with along controls
clean_adj_norm_medians %>%
  filter(!(Med_B2A_FSC<1.5 & Type == "Experimental")) %>%  #filter out outliers (likely resulting from contamination) as defined by Fluor <1.5
  ggplot(aes(generation, Med_B2A_FSC, color= sample)) +
  geom_line(aes(linetype = Type), size = 2.0) +
  scale_linetype_manual(values = c("dashed", "dashed", "dashed", "solid")) +
  scale_color_manual(values = c(
    "black", "black", "black", #controls
    "gray","gray","gray","gray","gray", #wildtype, 5, gray
    "#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59","#D7B02F","#CAA426","#D9BB59", #LTR,8,gold
    "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", #ARS, 7, softer salmon repeats
    "#6699cc","#6699cc","#6699cc","#6699cc","#6699cc","#6699cc","#6699cc","#6699cc" #LTR,8,
  )) +
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
#ggsave("MedNormFluo_FacetPlots_NoOutliers_010722.png")
#ggsave("MedNormFluro_v2_010722.png")
ggsave("MedNormFluor_051022.png")

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


# Make ridgeplots for each population

### on HPC: For Loop - for each sample, subset it and write a sc_distributions_SampleName_allTimepoints.csv
#for(pop in unique(sc_distr_alltimepoints$sample)) {
#  print(pop)
#  sc_distr_alltimepoints %>%
#  filter(sample == pop) %>%
#  write_csv(paste0("sc_distributions_",pop,"_all_timepoints.csv"))
#}

####################################################
# STEP 9:  Quantify CNV dynamics (Lauer et al. 2018)
# Author: Julie
  # 1) First, calculate Tup, the generation at which CNVs are initially detected, (Lang et al. 2011 and Lauer et al. 2018)
  # To do that, calculate the false positive rate for CNV detection (threshold), which I will define as the median frequency of 1 copy control cells appearing in the two_copy_or_more gate and gate across generations 8-260 plus the interquartile range (IQR) if the distribution of one copy controls appearing in the CNV gate as NOT normal. If the distribution is normal, then use the mean plus one standard deviation like in Lauer et al. 2018. Like in Lauer et al. 2018, samples surpassing this threshold is considered to contain CNVs.

#CNV False Positive Rate is defined by the frequency of the 1 copy control strain appearing in the CNV gate which is called the Two_or_more copy gate.
CNV_false_pos_df = #freq %>%
  freq_and_counts %>%
  filter(Count>70000) %>%
  anti_join(fails) %>%
  filter(Type == "1_copy_ctrl") %>%
  filter(Gate %in% c("two_or_more_copy")) %>%
  select(Type, Strain, Description, generation, Gate, Frequency, Count)

# draw a histogram, see if distribution is normal by eye
hist(CNV_false_pos_df$Frequency) #looks normal but could be left skewed
  abline(v = median(CNV_false_pos_df$Frequency),col = "red",lwd = 1.5)
  abline(v = mean(CNV_false_pos_df$Frequency), col = "blue", lwd = 1.5)

#Test for normality
shapiro.test(CNV_false_pos_df$Frequency) #null hypothesis is that the distribution is normal. if p <0.05, then it rejects the null hypothesis and so the distribution is NOT normal.
    #W = 0.95715, p-value = 0.4886
# Null hypothesis cannot be rejected. Our distribution is normal. Use the mean + 1SD as the threshold value.
# mean = 4.29, sd = 2.32
thres_mean = mean(CNV_false_pos_df$Frequency) + sd(CNV_false_pos_df$Frequency) #6.615341

#Determine Tup for each population, the generation when CNVs first appear above the threshold frequency.
Tup_per_pop = #freq %>%
  freq_and_counts %>%
  filter(Count>70000) %>%
  anti_join(fails) %>%
  filter(Type == "Experimental", Gate == "two_or_more_copy", Frequency >= thres_mean) %>%
  select(Type, Strain, Description, sample, generation, Gate, Frequency) %>% #View()
  group_by(sample) %>%
  slice(which.min(generation))

#Tup_per_pop %>% write_csv(file = "01_02_04_v2_fw_Tup_per_pop.csv")
#Tup_per_pop = read_csv(file = "01_02_04_v2_fw_Tup_per_pop.csv")

ggplot(Tup_per_pop, aes(reorder(Description, -generation),generation, fill = Description)) +
  geom_boxplot(outlier.shape = NA) +
  xlab("Genotype") +
  scale_fill_manual(values=c("#e26d5c", "#DEBD52", "#6699cc", "gray"))+
  ylab("Generation of first CNV appearance") +
  scale_x_discrete(labels=c("Wildtype architecture","LTR KO","ARS KO","LTR and ARS KO"))+
  scale_y_continuous(breaks=c(8,20, 30, 40, 50, 60, max(Tup_per_pop$generation)))+
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(family="Arial", size = 16, color = "black"), #edit x-tick labels
        axis.text.y = element_text(family="Arial", size = 20, color = "black"),
        text = element_text(size=18))+
  geom_jitter(size = 2, alpha = 0.9, color = c(
    rep("black",5), #wildtype, 5, gray
    "#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59","#D7B02F","#CAA426","#D9BB59", #LTR and ARS gold
    "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", "#e26d5c", #ARS, 7, softer salmon repeats
  rep("black", 8) #LTR,8, #blue
))
#ggsave("01_02_04_v2_fw_Tup_boxplot_042022.png")
#ggsave("01_02_04_v2_fw_Tup_boxplot_051022.png")

# ANOVA to test for significance
# One Way ANOVA because there is only 1 independent variable, genotype.
# Anova Tut: https://www.scribbr.com/statistics/anova-in-r/
Tup_anova = aov(generation~Description, data = Tup_per_pop)
summary(Tup_anova)
#            Df  Sum Sq Mean Sq F value  Pr(>F)
#Description  3   4483  1494.3   6.772 0.00181 **
#Residuals   24   5296   220.7
# Conclusion: There IS a significant difference in the means. Genotype has has significant effect on Tup.

#Calculate Sup
# "Sup is the rate of increase in CNV abundance during the initial expansion of the CNV subpopulation" Lauer et al. 2018
# S1 Text. Calculation of CNV dynamics parameters.
# Graphic representation of linear fit (and corresponding R2 values)
# during initial population expansion of CNV alleles.
# Slope of the linear fit corresponds to the dynamics parameter Sup shown in
# Table 1 and was calculated for the original evolution experiment and the
# barcode experiment. Data and code used to generate these figures can be
# accessed in OSF: https://osf.io/fxhze/. CNV, copy number variant.

#Calculate natural log proportion of each population with CNV relative to that without CNV
ln_table = freq_and_counts %>%
  filter(Count>70000) %>%
  filter(Gate %in% c("two_or_more_copy"), Type == "Experimental") %>%
  anti_join(fails)  %>% #remove contaminated and outliers informed by population ridgeplots (above) and fluor lineplots (below)
  group_by(sample, generation) %>%
  mutate(prop_CNV = sum(Frequency),
         prop_NoCNV = 100-prop_CNV,
         CNV_NoCNV = prop_CNV/prop_NoCNV,
         logECNV_NoCNV = log(CNV_NoCNV)) #log() function is natural logarithm in R (even though  log() commonly thought as base10 )

# JULIE: A function to calculate Sup, Explained Variance, make graphs, ggsave graphs
# then, use map() to apply function to all 28 populations - I have 28
# as part of the function, include a variable that changes number of fit points (window width)
# I recycled some code from Lauer et al. 2018

pop_list = unique(ln_table$sample) %>% sort()
gens = unique(ln_table$generation)

equation = function(x) {
  lm_coef <- list(a = round(coef(x)[1], digits = 2),
                  b = round(summary(x)[4]$coefficients[2], digits = 4),
                  r2 = round(summary(x)$r.squared, digits = 2));
  lm_eq <- substitute(slope == b~~~~italic(R)^2~"="~r2,lm_coef)
  as.character(as.expression(lm_eq));
}

#function to apply for loop to each of 28 populations, used some code from Lauer et al. 2018
sliding_fit = function(num_fitpoints, population){
  timepoints = nrow(subset(ln_table, sample %in% c(population)))
  rounds = timepoints- (timepoints/num_fitpoints) - 1
  m <- matrix(ncol = 6, nrow = rounds) #nrow = number of iterations. number of iterations depend on the number of generations and the number of fitpoints. max num of generations = 24. minimum num of fitpoints is 2. therefore nrow max is 23.
  colnames(m) <- c("start", "end", "gen_start", "gen_end", "slope", "rsquared")
  start = 1
  end = num_fitpoints
    for (i in 1:rounds){
      print(i)
      pop_data <- subset(ln_table, sample %in% c(population))
      fit_points_df <- subset(ln_table, sample %in% c(population) & generation >= gens[start] & generation <= gens[end])
      if (is.na(gens[end]) == TRUE ){
        break
      }
      fit <- lm(logECNV_NoCNV ~ generation, fit_points_df) #linear model, lm(y~x, by the data)
      print(summary(fit))

      ggplot(pop_data, aes(x=generation,y=(as.numeric(logECNV_NoCNV)), colour=sample)) +
      geom_point() +
      geom_smooth(data=fit_points_df, method=lm, show.legend=FALSE) +
      scale_y_continuous(expand = c(0, 0), 'ln(Prop. CNV/Prop. non-CNV)', limits = c(min(pop_data$logECNV_NoCNV)-1, max(pop_data$logECNV_NoCNV)+1)) +
      annotate("text", x = 200, y = min(pop_data$logECNV_NoCNV)-0.5, label = equation(fit), parse = TRUE) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 5), "Generations", limits=c(0,260)) +
      theme_classic() +
      scale_color_manual(values = c('black')) +
      guides(colour = guide_legend(override.aes = list(size=2))) +
      theme(legend.position = c(.15,.95), plot.title = element_text(size=14, hjust = 0.5), legend.title = element_blank(), axis.title.y = element_text(face="bold", size=12), axis.text.y = element_text(size=12), axis.title.x = element_text(face="bold", size=12), axis.text.x = element_text(size=12))

      ggsave(paste0(population,"_Sup_g",gens[start],"-",gens[end],"_",num_fitpoints,"pts.png"), width = 8, height = 5)

      print(paste0("From timepoints ",start," to ", end, ", generations ", gens[start], " to ", gens[end],", slope was ", as.numeric(coef(fit)[2]) %>% round(4)," and rsquared was ", as.numeric(summary(fit)[8]) %>% round(2) ))  #populate a data frame? five columns: start, end, gen_start, gen_end,  Rsq.

        m[i,1] <- start
        m[i, 2]<- end
        m[i, 3] <- gens[start]
        m[i,4] <- gens[end]
        m[i, 5] <- as.numeric(coef(fit)[2]) %>% round(4)
        m[i,6] <- as.numeric(summary(fit)[8]) %>% round(2)

        start = start + 1
        end = end + 1
        }
    m = m %>% na.omit() #remove NAs
    write_csv(as.data.frame(m), paste0(population,"_fits","_",num_fitpoints,"pts.csv"))
    return(m)
}

#call the function and assign outputted matrix to a variable
result <- sliding_fit(4, pop_list[27])
assign(paste0(pop_list[27],"_fits","_",4,"pts"), result) #Use assign() to rename and save to R environment

#call the function for all pops in the pop_list using map()
map(.x = pop_list[1:28], ~sliding_fit(4, .x)) #tilde inside functions (https://stackoverflow.com/questions/70665707/what-is-the-meaning-of-and-inside-the-function-map)

summary(fit)$coef[[4]] #fourth coefficient is the standard error of the linear model slope


####### MY PALLETTE

#Gold Metallic (6)
#DEBD52
#DBB741
#D7B02F
#CAA426
#D9BB59
#B89523

#NEW PALLETE 4/20/22
#WT = gray  "gray","gray","gray","gray","gray",
#LTR and ARS = GOLD metalic "#DEBD52","#DBB741","#D7B02F","#CAA426","#D9BB59","#D7B02F","#CAA426","#D9BB59", #LTR,8,gold
#ARS = SALMON  "#e26d5c"
#LTR = BABY BLUE "#6699cc"
