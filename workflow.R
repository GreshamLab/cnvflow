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

# Set working directory and get list of subdirectories
#setwd('../FCS files/') #Grace's working directory
#setwd('/Volumes/GoogleDrive/My Drive/Gresham Lab_Papers/2021/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files') #David's working directory
#setwd('G:/.shortcut-targets-by-id/1Bioj1YP_I7P8tqgmg4Zbt4EAfhb7J-0w/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files') #Titir's working directory
#setwd("/Volumes/GoogleDrive/My Drive/greshamlab/Molecular Determinants of CNV Evolution Dynamics/Summer 2021 Group LTEE/FCS files") #Julie's WD
setwd("~/Documents/Summer_LTEE_2021_FCS_files")
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
map(folders[23], make_exp_details, samplesheet = "EE_GAP1_ArchMuts_2021.csv")
#STEP 2: Read in all files in a directory and rename the channels.
#Results in one timepoint's gating set containing all .fcs files, associated experiment details, and marker details
#Author: Julie
folders = list.dirs()[-1]
exp_details_path = list.files(path = paste0(folders[1]), pattern = "_experiment_details.csv", full.names = T)
#exp_details_path = list.files(path = folder08, pattern = "_experiment_details.csv", full.names = T)
timept_gating_set <- cyto_setup(path = paste0(folders[1]), restrict=TRUE, select="fcs", details=F) #edit Markers on Viewer pane, Save & Close

#use pData to annotate the experiment details file associated with the gating set
experiment_details <- read_csv(exp_details_path) #import experiment-details.csv
flowWorkspace::pData(timept_gating_set)$name<-experiment_details$name
flowWorkspace::pData(timept_gating_set)$sample<-experiment_details$sample
flowWorkspace::pData(timept_gating_set)$`Outflow Well`<-experiment_details$`Outflow well`
flowWorkspace::pData(timept_gating_set)$Media<-experiment_details$Media
flowWorkspace::pData(timept_gating_set)$Strain<-experiment_details$Strain
flowWorkspace::pData(timept_gating_set)$Type<-experiment_details$Type
flowWorkspace::pData(timept_gating_set)$Description<-experiment_details$Description
flowWorkspace::pData(timept_gating_set)$generation<-experiment_details$generation

#file.rename(dir(pattern = "Experiment-Markers.csv"),"EE_GAP1_ArchMuts_2021-Experiment-Markers.csv") #rename the experiment-markers.csv file. Need to do once.

#STEP 3:  Perform gating on gating set
#Gate for 1) Cells, 2) Singlets, 3) CNVS
#Results in a gating file, and gates applied to all samples in the gating set.
#Author: Titir

#Log transform the data
# looks useful if I want to choose different transformation: https://dillonhammill.github.io/CytoExploreR/articles/CytoExploreR-Transformations.html
#timept_transformed <- cyto_transformer_log(timept_gating_set,
#                      channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A")) #returns it as a list
#transformed_timepoint_gating_set <- cyto_transform(timept_gating_set,
#                      trans = timept_transformed) #applies the the transformation and returns it as a gatingSet
#timept_transformed <- cyto_transformer_logicle(timept_gating_set,
#                                           channels = c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "B2-A"),
#                                           widthBasis = -10) #returns it as a list
#transformed_timepoint_gating_set <- cyto_transform(timept_gating_set,
#                                     trans = timept_transformed) #applies the the transformation and returns it as a gatingSet
#quickly check the transformation by plotting the data
cyto_plot_explore(transformed_timepoint_gating_set,
                  channels_x = "FSC-A",
                  channels_y = "GFP",
                  axes_limits = "data")

##Gating using the entire timepoint dataset.
#First we gate for the cells
cyto_gate_draw(transformed_timepoint_gating_set,
  #transformed_logicle_timept,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"),
               axes_limits = "data",
               gatingTemplate = "cytek_gating_test.csv",
)

#Then we define the singlets based on forward scatter height and width
cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Cells",
               alias = "Single_cells",
               channels = c("FSC-A","FSC-H"),
               axes_limits = "data",
               gatingTemplate = "cytek_gating_test.csv"
)

#Gating for CNVs using the 0,1 and 2 copy controls:
zero_copy <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(30,61)] #DGY1

one_copy <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(1,32)] #DGY500

two_copy <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(31,62)] #DGY1315

cyto_gate_draw(transformed_timepoint_gating_set,
               parent = "Single_cells", #first color
               alias = c("zero_copy", "one_copy", "two_copy","multi_copy"), #defines gate names
               channels = c("FSC-A","B2-A"),
               axes_limits = "data",
               select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
               #gatingTemplate = "cytek_gating_test.csv",
               overlay = c(zero_copy, one_copy, two_copy),
               point_col = c("black", "green", "red", "blue")
)

#STEP 4:  Generate statistics tables
#Results in a .csv file in tidy format that includes all metadata and specifies proportion of cells with 0, 1, 2, 3+ copies, a .csv file with overall median GFP and FSC-A, and a third .csv file with gate-wise median GFP and FSC-A.
#Author: Titir

stats_freq_01 <- cyto_stats_compute(transformed_timepoint_gating_set01,
                                  parent = "Single_cells",
                                  alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
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
                                              alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                              channels = c("FSC-A", "B2-A"),
                                              stat="median",
                                              save_as = "stats_median_gatewise_01.csv")

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

  my_path <- folder_name #gets relative path name for folder to be analyzed

  prefix <- folder_name %>% str_extract("([0-9])+_EE_GAP1_ArchMuts_2021") #extracts the time point number from folder name

  my_expt_details_path <- paste0(my_path,"/",prefix,"_experiment_details.csv") #gets experiment details .csv from correct directory
  #my_expt_details_path <- list.files(path = paste0(my_path), pattern = "_experiment_details.csv", full.names = T)
  #1. read in files and make a gating set
  print(my_path)
  timepoint_gating_set <- cyto_setup(path=my_path, select="fcs", details=F, markers = F)

  #2. read in experiment details for that gating set
  my_experiment_details <- read_csv(my_expt_details_path, show_col_types = F) #import experiment-details.csv
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
  #               gatingTemplate = "cytek_gating_test.csv",
  #)

  #Then we define the singlets based on forward scatter height and width
  #cyto_gate_draw(transformed_timepoint_gating_set,
  #               parent = "Cells",
  #               alias = "Single_cells",
  #               channels = c("FSC-A","FSC-H"),
  #               axes_limits = "data",
  #               gatingTemplate = "cytek_gating_test.csv"
  #)
#  zero_copy <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(30,61)] #DGY1
#  one_copy <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(1,32)] #DGY500
#  two_copy <- cyto_extract(transformed_timepoint_gating_set, "Single_cells")[c(31,62)] #DGY1315
#  cyto_gate_draw(transformed_timepoint_gating_set,
#                 parent = "Single_cells", #first color
#                 alias = c("zero_copy", "one_copy", "two_copy","multi_copy"), #defines gate names
#                 channels = c("FSC-A","B2-A"),
#                 axes_limits = "data",
#                 select = list(Strain = c("DGY1","DGY500","DGY1315")),  #control strains
#                 gatingTemplate = "cytek_gating_test.csv",
#                 overlay = c(zero_copy, one_copy, two_copy),
#                 point_col = c("black", "green", "red", "blue")
#  )

  #5. apply gating-template.csv to transformed gating set
  cyto_gatingTemplate_apply(transformed_timepoint_gating_set, gatingTemplate= gating_template)

  #6. write stats: freq file for % of cells inside each gate, median FSC and GFP for each population, median FSC and GFP for each gated population
  #Titir
  stats_freq <- cyto_stats_compute(transformed_timepoint_gating_set,
                                      parent = c("Single_cells"),
                                     alias = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                      stat="freq",
                                      save_as = paste0("01_11_GFPlogicle_FSClog_v3_stats_freq_",prefix,".csv") #writes to working directory
                                      )
 stats_median_overall <- cyto_stats_compute(transformed_timepoint_gating_set,
                                     parent = c("Single_cells"),
                                     alias  = c("Single_cells"),
                                     channels = c("FSC-A", "B2-A"),
                                     stat="median",
                                   save_as = paste0("01_11_GFPlogicle_FSClog_v3_stats_median_overall_", prefix,".csv"))

  stats_cell_number <- cyto_stats_compute(transformed_timepoint_gating_set,
                                             parent = c("Single_cells"),
                                             alias  = c("Single_cells"),
                                             channels = c("FSC-A", "B2-A"),
                                             stat="count",
                                             save_as = paste0("01_11_GFPlogicle_FSClog_v3_stats_cell_number_", prefix,".csv"))

 stats_median_gatewise <- cyto_stats_compute(transformed_timepoint_gating_set,
                                              parent = c("Single_cells"),
                                              alias  = c("zero_copy", "one_copy", "two_copy", "multi_copy"),
                                              channels = c("FSC-A", "B2-A"),
                                              stat="median",
                                             save_as = paste0("01_11_GFPlogicle_FSClog_v3_stats_median_gatewise_", prefix,".csv"))
}

#STEP 6:  Apply function from STEP 5 to all subdirectories
#Uses map from purr() to apply function from step 5 to all directories
#Author: Julie

#map(folders[-1], analyze_all_exp, my_markers, gating_template = "cytek_gating.csv")
try(map(folders[18:length(folders)],analyze_all_exp, my_markers, gating_template = "cytek_gating_test.csv"))
try(map(folders[23],analyze_all_exp, my_markers, gating_template = "cytek_gating_test.csv"))
#STEP 7:  Combine stats_freq.csv and stats_median.csv files into a single dataframe
#Pull in all stats_* files from directories and assemble into a single dataframe
#Author: Julie

list.files(path = ".", pattern = "01_11_GFPlogicle_FSClog_v3_stats_freq") %>%
  read_csv() %>%
  write_csv(file = "01_11_GFPlogicle_FSClog_v3_stats_freq_all_timepoints.csv")

list.files(path = ".", pattern = "01_11_GFPlogicle_FSClog_v3_stats_median_overall") %>%
  read_csv() %>%
  write_csv(file = "01_11_GFPlogicle_FSClog_v3_stats_median_overall_all_timepoints.csv")

list.files(path = ".", pattern = "01_11_GFPlogicle_FSClog_v3_stats_median_gatewise") %>%
  read_csv() %>%
  write_csv(file = "01_11_GFPlogicle_FSClog_v3_stats_median_gatewise_all_timepoints.csv")

list.files(path = ".", pattern = "01_11_GFPlogicle_FSClog_v3_stats_cell_number") %>%
  read_csv() %>%
  write_csv(file = "01_11_GFPlogicle_FSClog_v3_stats_cell_number_all_timepoints.csv")

#STEP 8: Plot time series & assess gates
#Determine whether =>83% of controls are in the correct gate
#Author: Grace & Julie

# read in frequency csv, median csvs for all timepoints
freq = read_csv("01_11_GFPlogicle_FSClog_v3_stats_freq_all_timepoints.csv") %>% rename(Gate = Population)
medians = read_csv("01_11_GFPlogicle_FSClog_v3_stats_median_overall_all_timepoints.csv")
medians_bygate = read_csv("01_11_GFPlogicle_FSClog_v3_stats_median_gatewise_all_timepoints.csv")
cell_numbers = read_csv("01_11_GFPlogicle_FSClog_v3_stats_cell_number_all_timepoints.csv")

# add cell number column to freq table
freq = left_join(freq, cell_numbers) %>%
  select(-Marker)

# exclude any well/timepoint with less than 70,000 single cells
fails = freq %>%
  filter(Count>70000) %>%
# check controls are in their proper gates
  filter(str_detect(Description, "control")) %>%
  select(Description, Strain, generation, Gate, Frequency, name, Count) %>%
  mutate(flag = case_when(Strain == "DGY1" & Gate == "zero_copy" & Frequency >= 83 ~ "pass",
                          Strain == "DGY1" & Gate == "zero_copy" & Frequency < 83 ~ "fail",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency >= 83 ~ "pass",
                          Strain == "DGY500" & Gate == "one_copy" & Frequency < 83 ~ "fail",
                          Strain == "DGY500" & Gate == "zero_copy" & Frequency >= 15 ~ "fail",
                          Strain == "DGY1315" & Gate == "two_copy" & Frequency >= 83 ~ "pass",
                          Strain == "DGY1315" & Gate == "two_copy" & Frequency < 83 ~ "fail",
                          Strain == "DGY1315" & Gate == "zero_copy" & Frequency >= 15 ~ "fail",
                          Strain == "DGY500" & Gate == "two_copy" & Frequency >=15 ~ "fail"
                          ))%>%
  filter(flag == "fail") %>%
  arrange(Description)
  View(fails)
  fails %>% write_csv("01_11_GFPlogicle_FSClog_v3_83_fail.csv")

# plot controls over time
freq %>%
filter(Count>70000) %>%
  filter(str_detect(Description, "control")) %>%
  filter(Type == "0_copy_ctrl", generation >200) %>%
  select(Type, Strain, Description, generation, Gate, Frequency, Count) %>% View()
  #anti_join(fails) %>% #exclude the contaminated controls timepoints (the failed timepoints)
  ggplot(aes(generation, Frequency, color = Gate)) +
  geom_line() +
  facet_wrap(~Description) +
  ylab("% of cells in gate") +
  theme_minimal() +
  theme(text = element_text(size=16))

# plot proportion of population in each gate over time for all experimental
#JULIE: edit to exclude the contaminated control timepoints
plot_list = list()
i=1
for(exp in unique(freq$Description)) {
  #print(plot_dist(obs))
  plot_list[[i]] = freq %>%
    filter(Count>70000) %>%
    filter(Description==exp) %>%
    ggplot(aes(generation, Frequency, color = Gate)) +
    geom_line() +
    facet_wrap(~sample) +
    ylab("% of cells in gate") +
    theme_minimal()+
    theme(text = element_text(size=18))
  i = i+1
}
names(plot_list) = unique(freq$Description)
plot_list$`GAP1 WT architecture` # change index to view replicates for different genetic backgrounds
plot_list$`GAP1 ARS KO`
plot_list$`GAP1 LTR KO`
plot_list$`GAP1 LTR + ARS KO`

# plot proportion of the population with a CNV over time
#Julie: dont need controls
freq %>%
  filter(Count>70000) %>%
  #anti_join(fails) %>%
  group_by(sample, generation) %>% #View()
  filter(Gate %in% c("two_copy", "multi_copy"), Type == "Experimental") %>%
  group_by(sample, generation) %>%
  mutate(prop_CNV = sum(Frequency)) %>% #View()
  select(sample, generation, Description, prop_CNV) %>%
  distinct() %>%
  ggplot(aes(generation, prop_CNV, color = sample)) +
  geom_line() +
  facet_wrap(~Description) +
  theme_minimal() +
  ylab("Proportion of the population with GAP1 CNV") +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme(text = element_text(size=16), legend.position = "none")

#Quantify CNV dynamics
  # 1) First, calculate Tup, the generation at which CNVs are initially detected, (Lang et al. 2011 and Lauer et al. 2018)
  # To do that, calculate the false positive rate for CNV detection (threshold), which I will define as the median frequency of one-copy control cells appearing in the two copy and multicopy gate across generations 8-260.
  # In Lauer et al. 2018, the false positive for CNV detection is defined as the average plus one standard deviation.
  # Since my distribution is not normal, I will use the median plus IQR instead as my false positive threshold.
  # Like in Lauer et al. 2018, samples surpassing this

CNV_false_pos_freq = freq %>%
  filter(Count>70000) %>%
  anti_join(fails) %>%
  filter(generation != 174, Type == "1_copy_ctrl") %>%
  filter(Gate %in% c("two_copy", "multi_copy")) %>%
  group_by(generation) %>%
  mutate(CNV_freq = as.numeric(sum(Frequency))) %>%
  filter(Gate == "multi_copy")
 mean(CNV_false_pos_freq$CNV_freq)
 sd(CNV_false_pos_freq$CNV_freq)

 hist(CNV_false_pos_freq$CNV_freq) #not normal
  abline(v = median(CNV_false_pos_freq$CNV_freq),col = "red",lwd = 1.5)
  abline(v = mean(CNV_false_pos_freq$CNV_freq), col = "blue", lwd = 1.5)
 m = median(CNV_false_pos_freq$CNV_freq) #3.96294
 iqr = IQR(CNV_false_pos_freq$CNV_freq) #3.159226
 threshold = m + iqr #threshold from median + IQR #7.122166
 mean(CNV_false_pos_freq$CNV_freq) #4.113866
 sd(CNV_false_pos_freq$CNV_freq) #2.72
 mean(CNV_false_pos_freq$CNV_freq) + sd(CNV_false_pos_freq$CNV_freq) #threshold from mean + 1SD #6.81973

#frequency of experimental sample cells in two_copy and multi_copy gates over generation
 # sooo this is the same as the propCNV plot
freq %>%
   filter(Count>70000) %>%
   filter(generation != 174, Type == "Experimental") %>%
   #anti_join(fails) %>% #exclude the contaminated controls timepoints (the failed timepoints)
   filter(Gate %in% c("two_copy", "multi_copy")) %>%
  group_by(sample, generation) %>%
  mutate(CNV_freq = as.numeric(sum(Frequency))) %>%
  filter(Gate == "multi_copy") %>%
   ggplot(aes(generation, CNV_freq, color = sample)) +
   geom_line() +
   facet_wrap(~Description) +
   ylab("% of cells in gate") +
   theme_minimal() +
   theme(text = element_text(size=16))

#same as above but the median CNV_frequency of the populations


# plot ridgeplots (histograms):
# normalized fluorescence histograms of controls at each generations like in Lauer et al. Fig 2A.
# aka not looking at the median GFP values per population per generation.
  # We want to see the entire distribution of GFP cells per population per generation.


# plot median proportion of the populations with a CNV over time
  freq %>%
    filter(Count>70000) %>%
    group_by(sample, generation) %>%
    #filter(generation != 174, sample != "gap1_4") %>%
    filter(Gate %in% c("two_copy", "multi_copy"), Type == "Experimental") %>%
    group_by(sample, generation) %>%
    mutate(prop_CNV = sum(Frequency)) %>%
    arrange(generation, Description) %>%
    group_by(Description, generation) %>%
    mutate(median_propCNV = median(prop_CNV),
           mean_propCNV = mean(prop_CNV)) %>% #View()
      filter(is.na(generation)) %>% View()
  select(sample, generation, Description, median_propCNV) %>%
    distinct() %>% #View()
    ggplot(aes(generation, median_propCNV, color = Description)) +
    geom_line() +
    facet_wrap(~Description) +
    theme_minimal() +
    ylab("Median proportion of the population with GAP1 CNV") +
    scale_x_continuous(breaks=seq(0,250,50)) +
    theme(text = element_text(size=14), legend.position = "none")

#plot median GFP fluorescence normalized over median FSC-A over time for experimental
  #overlay 0,1,2 controls on same graph as experimental with gray lines
norm_medians = medians %>%
  pivot_wider(names_from = Marker, values_from = MedFl) %>%
  mutate(NormMedGFP = GFP/`FSC-A`) %>% #normalization
  arrange(Description) %>%
  slice(rep(1:56, each = 4)) %>% #repeat each control row 4 times
  mutate(Description = rep(c("GAP1 ARS KO", "GAP1 LTR + ARS KO", "GAP1 LTR KO","GAP1 WT architecture"), times=56)) %>% #Mutate Description of controls so we can graph them on experimental facet plots
  merge(norm_medians %>% filter(Type == "Experimental"), all = TRUE) %>% #merge back to experimental rows
  arrange(Description, generation) %>%
  filter(is.na(generation)) %>% View()

#graph experimental with along controls
ggplot(norm_medians, aes(generation, NormMedGFP, color= sample)) +
  geom_line(aes(linetype = Type)) +
  scale_linetype_manual(values = c("dashed", "dashed", "dashed", "solid")) +
  #scale_color_manual(values=c("gray", "gray", "gray", rep("orange", 5), rep("purple", 8),rep("blue", 7),rep("pink", 8) ))+
  #scale_color_brewer(palette = "Set1")+
  #scale_fill_brewer(value(-1, "black")+
  facet_wrap(~Description) +
  ylab("normalized median fluorescence") +
  scale_x_continuous(breaks=seq(0,250,50)) +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size=16))

# plot median of the populations' median normalized fluorescence over time
# aka collapse the populations - take the median across the populations
norm_medians %>%
  group_by(generation, Description) %>%
    filter(Type == "Experimental" ) %>%
    mutate(pops_median = median(NormMedGFP) ) %>%
  filter(is.na(generation)) %>% View()




