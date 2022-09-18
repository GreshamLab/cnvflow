###### Plot normalized median mCitrine fluorescence over time with controls on same plot
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
  left_join(cell_numbers) %>%

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
  anti_join(adj_norm_medians %>% filter)
anti_join(adj_norm_medians %>% filter(generation == 231 & Type == "0_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 182 & Type == "1_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 203 & Type == "1_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 252 & Type == "1_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 260 & Type == "1_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 79 & Type == "2_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 95 & Type == "2_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 108 & Type == "2_copy_ctrl")) %>%
  anti_join(adj_norm_medians %>% filter(generation == 116 & Type == "2_copy_ctrl")) %>%

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


### # Graph single plots of median normalized fluorescence for every population (sample)
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

