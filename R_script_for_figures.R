##### Variations in Atlantic water influx and sea-ice cover drive taxonomic and 
##### functional shifts in Arctic marine bacterial communities

### R script to process data tables and reproduce figures ###

### It is essential that the first two sections of the script are run
### 'LOADING PACKAGES AND DATA' and 'PREFILTERING ASV DATA AND ASSESSING 
### DISTRIBUTION DYNAMICS'
### These sections include necessary processing of the data and generate
### dataframes that are then used throughout.
### Once these first sections are run, all other individual Figure sections can
### be run independently

### Any errors or inconsistencies, feel free to contact tpriest@mpi-bremen.de ###

### Set working directory!!!

setwd("C:/Users/Taylor/ownCloud/PhD/Arctic Bacterioplankton/Fram Strait/RAS_metagenomes/EGC/ASV/FRAM_RAS_EGC_data_for_processing")


############################################################################################
### LOADING PACKAGES AND DATA ###
############################################################################################

# Loading libraries and data tables
library(mixOmics)
library(phyloseq)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(factoextra)
library(GGally)
library(textshape)
library(ggord)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(gtools)
library(amap)
library(ape)
library(ALDEx2)
library(rgr)
library(reshape2)

# Load ASVs
ASV <- read.table(
  "FRAM_RAS_EGC_ASV_raw.txt",
  h = T, sep = "\t",
  check.names=F,
  row.names=1)

# Load taxonomy
TAX <- read.table(
  "FRAM_RAS_EGC_ASV_taxa.txt",
  h = T, 
  sep = "\t",
  check.names=F, 
  row.names=1)
TAX <- as.matrix(TAX)

# Load Metadata
# Assign ice/PW categories 
ENV <- read.table(
  "FRAM_RAS_EGC_ASV_meta.txt", header=T) %>%
  mutate(Mooring_position=case_when(Mooring_full %in% c(
    "EGC-3","EGC-4")~"MIZ",
    TRUE~"core-EGC")) %>%
  mutate(Ice_cover_category=case_when(
    Ice_cover > 70 ~"high",
    Ice_cover > 40 & Ice_cover < 70 ~"medium",
    Ice_cover < 40 ~"low")) %>%
  mutate(Ice_cover_past_category=case_when(
    Ice_cover_past > 70 ~"high",
    Ice_cover_past > 40 & Ice_cover_past < 70 ~"medium",
    Ice_cover_past < 40 ~"low")) %>% 
  mutate(Mooring_name = paste(Mooring_position,Year,sep="_"))

############################################################################################
### PREFILTERING ASV DATA AND ASSESSING DISTRIBUTION DYNAMICS ###
############################################################################################

### CREATING PHYLOSEQ OBJECT, PREFILTERING AND CONVERTING TO RELATIVE ABUNDANCE

# Reformat
asv = otu_table(ASV, taxa_are_rows=T)
tax = tax_table(TAX)
rownames(tax) <- rownames(asv)
env = ENV %>%
  column_to_rownames("RAS_id")

# Subset: ASVs with >3 counts in >3 samples
pseq.abs <- phyloseq(
  otu_table(asv, taxa_are_rows = T), 
  sample_data(env), 
  tax_table(tax)) %>%
  phyloseq::filter_taxa(function(x) sum(x >= 3) >= 3, TRUE)

# Fix tax-IDs
colnames(pseq.abs@tax_table) <- c(
  "Kingdom","Phylum","Class","Order",
  "Family","Genus","Species")

# Rel. abundances
pseq.rel = transform_sample_counts(
  pseq.abs, 
  function(x) x / sum(x) * 100) 

# Extract filtered ASV dataset with absolute and relative matrices

# Absolute counts
asv_abs <- as(otu_table(pseq.abs), "matrix") %>% as.data.frame() 

# Relative abundances
asv_rel <- as(otu_table(pseq.rel), "matrix") %>% as.data.frame() %>% mutate_all(round, 3)

# Taxonomy
asv_tax <- as(tax_table(pseq.abs), "matrix") %>%  as.data.frame()

# Save filtered datasets
write.table(asv_abs, file="FRAM_RAS_EGC_ASV_absolute_filtered.txt", sep="\t")
write.table(asv_rel, file="FRAM_RAS_EGC_ASV_relative_filtered.txt", sep="\t")
write.table(asv_tax, file="FRAM_RAS_EGC_ASV_taxonomy_filtered.txt", sep="\t")
write.table(ENV, file="FRAM_RAS_EGC_ASV_meta_refined.txt", sep="\t")

### ASSESSING DISTRIBUTION DYNAMICS OF ASVS

# Calculate frequency of detection across samples
asv_freq_detection = asv_abs %>% 
  replace(., . > 1, 1) %>% 
  mutate(Num_samples_present = rowSums(.[1:84])) %>% 
  as.data.frame(.) %>% 
  tibble::rownames_to_column(., var="ASV_name") %>%
  subset(., select=c(ASV_name,Num_samples_present))

# Calculate whether ASV present at both moorings or only one
asv_mooring_distribution = asv_abs %>% 
  replace(., . > 1, 1) %>% 
  tibble::rownames_to_column(., var="ASV_name") %>% 
  reshape2::melt(., id.vars="ASV_name", variable.name="RAS_id", value.name="Rel_abund") %>%
  left_join(ENV, by="RAS_id") %>% 
  subset(., select=c(ASV_name,Mooring_position,Rel_abund)) %>%
  aggregate(Rel_abund~ASV_name+Mooring_position, data=., FUN=sum) %>% 
  reshape2::dcast(., ASV_name~Mooring_position, value.var="Rel_abund") %>% 
  mutate(Mooring_presence = case_when(
    `core-EGC` > 0 & MIZ > 0 ~ "Shared", 
    `core-EGC` >= 1 & MIZ < 1 ~ "core-EGC", 
    `core-EGC` < 1 & MIZ > 0 ~ "MIZ")) %>% 
  subset(., select=c(ASV_name,Mooring_presence))

# Maximum relative abundance of each ASV
asv_max_rel_abund = asv_rel %>% 
  tibble::rownames_to_column(., var="ASV_name") %>%
  melt(., id.vars="ASV_name", variable.name="Sample", value.name="Rel_abund") %>% 
  aggregate(Rel_abund~ASV_name, data=., max) %>% 
  rename(., Max_rel_abund = Rel_abund)

# Min relative abundance of each ASV
asv_min_rel_abund = asv_rel %>% 
  tibble::rownames_to_column(., var="ASV_name") %>%
  melt(., id.vars="ASV_name", variable.name="Sample", value.name="Rel_abund") %>% 
  aggregate(Rel_abund~ASV_name, data=., min) %>% 
  rename(., Min_rel_abund = Rel_abund)

# Average relative abundance of each ASV
asv_avg_rel_abund = asv_rel %>% 
  tibble::rownames_to_column(., var="ASV_name") %>%
  melt(., id.vars="ASV_name", variable.name="Sample", value.name="Rel_abund") %>% 
  aggregate(Rel_abund~ASV_name, data=., mean) %>% 
  rename(., Avg_rel_abund = Rel_abund)

# Average relative abundance of each ASV at both moorings and classifiying mooring
# preference based on the higher of these values
asv_mooring_pref_and_avg_rel_abund_at_moorings = asv_rel %>% 
  tibble::rownames_to_column(., var="ASV_name") %>%
  melt(., id.vars="ASV_name", variable.name="Sample", value.name="Rel_abund") %>%
  rename(., RAS_id = Sample) %>%
  left_join(ENV, by="RAS_id") %>% 
  subset(., select=c(ASV_name,Mooring_position,Rel_abund)) %>% 
  aggregate(Rel_abund~ASV_name+Mooring_position, data=., mean) %>% 
  dcast(., ASV_name~Mooring_position, value.var="Rel_abund") %>% 
  mutate(Pref_region=case_when(
    `core-EGC` > MIZ ~ "core-EGC",
    MIZ > `core-EGC` ~ "MIZ")) %>% 
  as.data.frame(.) %>%
  rename(., `Avg_rel_abund_core-EGC` = `core-EGC`) %>%
  rename(., Avg_rel_abund_MIZ = MIZ)

# Max relative abundance of each ASV at both moorings
asv_max_rel_abund_at_moorings = asv_rel %>% 
  tibble::rownames_to_column(., var="ASV_name") %>%
  melt(., id.vars="ASV_name", variable.name="Sample", value.name="Rel_abund") %>%
  rename(., RAS_id = Sample) %>%
  left_join(ENV, by="RAS_id") %>% 
  subset(., select=c(ASV_name,Mooring_position,Rel_abund)) %>% 
  aggregate(Rel_abund~ASV_name+Mooring_position, data=., max) %>% 
  dcast(., ASV_name~Mooring_position, value.var="Rel_abund") %>% 
  as.data.frame(.) %>%
  rename(., `Max_rel_abund_core-EGC` = `core-EGC`) %>%
  rename(., Max_rel_abund_MIZ = MIZ)

### Create a master ASV dataframe containing all of the above-calculated 
### statistics along with taxonomy

tax_temp <- tibble::rownames_to_column(asv_tax, var="ASV_name")

asv_summary_temp = asv_freq_detection %>% 
  left_join(asv_max_rel_abund, by="ASV_name") %>%
  left_join(asv_min_rel_abund, by="ASV_name") %>%
  left_join(asv_avg_rel_abund, by="ASV_name") %>%
  left_join(asv_mooring_distribution, by="ASV_name") %>% 
  left_join(asv_mooring_pref_and_avg_rel_abund_at_moorings, by="ASV_name") %>% 
  left_join(asv_max_rel_abund_at_moorings, by="ASV_name") %>%
  left_join(tax_temp, by="ASV_name")

### Classify ASVs into resident, intermittent and transient based on frequency
### of detection

asv_summary_complete = asv_summary_temp %>% 
  mutate(Distribution_group = case_when(
    Num_samples_present > 76 ~ "Resident",
    Num_samples_present > 21 & Num_samples_present < 77 ~ "Intermittent",
    Num_samples_present < 22 ~ "Transient"
  ))

# Export dataframe
write.table(asv_summary_complete, file="FRAM_RAS_EGC_ASV_dynamics_summary_and_taxa_complete.txt", 
            sep="\t")


############################################################################################
### Figure 1 ###
############################################################################################

##########

### Figure 1. Geographical location of seafloor moorings and variation in 
### environmental conditions in MIZ (2016-2018) and core-EGC (2018-2020).

##########

# Part A and B for Figure 1 were generated using QGIS 

# Part C:

# Import data and subset
egc_asv_meta <- read.table(
  "FRAM_RAS_EGC_ASV_meta_refined.txt", header=T)

# Subset data and rename columns
egc_meta_subset = egc_asv_meta %>% 
  subset(., select=c(Mooring_name,Temperature,AW_proportion,Ice_cover)) %>% 
  melt(id.vars=c("Mooring_name"), 
       variable.name="variable", value.name="value") %>% 
  mutate(variable = case_when(
    (variable == "Temperature" ~ "Temperature (°C)"),
    (variable == "AW_proportion" ~ "AW proportion (%)"),
    (variable == "Ice_cover" ~ "Ice cover (%)")))

# Plot Figure 1c
egc_mooring_metadata_boxplot <- ggplot(egc_meta_subset) + 
  geom_boxplot(aes(x=Mooring_name,y=value,fill=Mooring_name), stat="boxplot") + 
  facet_wrap(.~variable, scales="free_y") + 
  labs(fill = "Mooring") + 
  theme_bw() + 
  scale_fill_manual(values=c("#db6d00","#924900","#009292","#004949")) + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.background.x = element_rect(fill = "white", colour = "black"))

# Export Figure 1c
pdf(file = "Figure_1_part_C.pdf", height=5, width=10)
egc_mooring_metadata_boxplot
dev.off()

# Figure 1a, 1b and 1c were then combined in Inkscape

############################################################################################
### Figure 2 ###
############################################################################################

##########

### Figure 2. Community structure across variations in water mass, ice cover and 
### daylight conditions.

##########

# Import ASV_relative_abundance_matrix
egc_asv_rel=read.table("FRAM_RAS_EGC_ASV_relative_filtered.txt", header = TRUE, 
                       sep = "\t",as.is=TRUE,check.names=F,row.names=1)

# Import sample metadata
egc_asv_meta=read.table("FRAM_RAS_EGC_ASV_meta_refined.txt", header = TRUE, 
                        sep = "\t",check.names=F,as.is=TRUE)

# Reformat data and calculate Bray-Curtis dissimilarity matrix
egc_asv_bray_dist = egc_asv_rel %>% 
  t(.) %>% 
  as.data.frame(.) %>% 
  tibble::rownames_to_column(., var="Sample") %>% 
  arrange(., Sample) %>% 
  textshape::column_to_rownames(., "Sample") %>% 
  as.matrix(.) %>% 
  vegdist(., distance="bray", na.rm=TRUE)

# Reformat metadata and standardize 
egc_asv_meta[is.na(egc_asv_meta)]<-0

egc_asv_meta_stand <- egc_asv_meta %>% 
  subset(., select=c(RAS_id,Temperature,Salinity,O2_concentration,AW_proportion,Daylight,
                     Ice_cover,Ice_edge_distance,Ice_cover_past)) %>% 
  arrange(., RAS_id) %>%
  textshape::column_to_rownames(., "RAS_id") %>%
  decostand(., method="standardize", MARGIN=2)

# Check for colinearity amonst metadata variables
ggpairs(egc_asv_meta_stand)

# Significant correlations >0.85 are considered highly colinear
# Temperature and salinity are highly colinear with AW_proportion
# Will remove temperature and salinity and keep AW proportion

egc_asv_meta_stand_filt <- subset(egc_asv_meta_stand, select=-c(Temperature,Salinity))

# Distance-based redundancy analysis with ordistep selection of variables
egc_asv_dbrda_null <- dbrda(egc_asv_bray_dist~1, data = egc_asv_meta_stand_filt,
                            sqrt.dist=TRUE, add=TRUE)
egc_asv_dbrda <- dbrda(egc_asv_bray_dist~., data = egc_asv_meta_stand_filt,
                       sqrt.dist=TRUE, add=TRUE)

egc_asv_dbrda_ordistep <- ordiR2step(egc_asv_dbrda_null, 
                                     scope = formula(egc_asv_dbrda), 
                                     direction="both",
                                     pstep=1000, r2scope=TRUE, permutations=1000)
egc_asv_dbrda_ordistep

# Test the significance of the model
anova.cca(egc_asv_dbrda_ordistep, step=1000)

# Select these variables and rerun the dbRDA
egc_asv_dbrda_significant <- dbrda(egc_asv_bray_dist~AW_proportion+Daylight+
                                     Ice_cover_past,
                                   data = egc_asv_meta_stand_filt)

# Check proportion of constrained variation
egc_asv_dbrda_significant

# Test model significance and axes significance 
anova.cca(egc_asv_dbrda_significant, step=1000)
anova.cca(egc_asv_dbrda_significant, by='axis', step=1000)
anova.cca(egc_asv_dbrda_significant, by='margin', step=1000, model="direct")
vif.cca(egc_asv_dbrda_significant)

# Calculate adjusted R-squared value of model
RsquareAdj(egc_asv_dbrda_significant)$adj.r.squared

# Extract sample scores
egc_asv_dbrda_sample_scores <- egc_asv_dbrda_significant$CCA$wa

# Convert output into ggplot object
temp_ord <- ggord(egc_asv_dbrda_significant)
temp_ord

# Extract information from plot and reformat data
egc_asv_dbrda_ord_meta <- temp_ord$plot_env$vecs
egc_asv_dbrda_ord_meta_ref = egc_asv_dbrda_ord_meta %>% 
  rename(xend = one) %>% 
  rename(yend = two) %>% 
  rename(Variable = lab) %>%
  mutate(xbeg = 0,
         ybeg = 0)

egc_asv_dbrda_ord_points <- temp_ord$plot_env$obs
egc_asv_dbrda_ord_points_and_meta = egc_asv_dbrda_ord_points %>% 
  rename(CAP1 = one) %>% 
  rename(CAP2 = two) %>% 
  rename(RAS_id = lab) %>% 
  left_join(egc_asv_meta, by = "RAS_id")
egc_asv_dbrda_ord_points_and_meta

### Plotting Figure 2

# Plot Figure 2a (AW proportion)
egc_asv_dissim_pca_aw_proportion = ggplot(egc_asv_dbrda_ord_points_and_meta, aes(CAP1, CAP2)) + 
  geom_point(aes(colour=AW_proportion*100, shape=Mooring_position), size=4) + 
  scale_colour_gradient2(low = "#D9D9D9", high="#08204F",
                         mid="#6BAED6", midpoint=50, 
                         guide="colourbar") +  
  labs(x = "CAP1 (16%)", y = "CAP2 (7%)", colour="AW proportion (%)") + 
  theme_bw() + 
  theme(legend.position="bottom", 
        legend.title = element_text(size = 14, colour="black"), 
        legend.text = element_text(size = 12, colour="black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.title.x = element_text(size = 14, colour = "black")) + 
  guides(shape=guide_legend(override.aes = list(size=5)))

# Plot Figure 2b (past ice cover)
egc_asv_dissim_pca_ice = ggplot(egc_asv_dbrda_ord_points_and_meta, aes(CAP1, CAP2)) + 
  geom_point(aes(colour=Ice_cover_past, shape=Mooring_position), size=4) + 
  scale_colour_gradient2(low = "#D9D9D9", high="#0f2310",
                         mid="#91c591", midpoint=50, 
                         guide="colourbar") + 
  labs(x = "CAP1 (16%)", y = "CAP2 (7%)", colour="Past ice cover (%)") + 
  theme_bw() + 
  theme(legend.position="bottom", 
        legend.title = element_text(size = 14, colour="black"), 
        legend.text = element_text(size = 12, colour="black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.title.x = element_text(size = 14, colour = "black")) + 
  guides(shape=guide_legend(override.aes = list(size=5)))

# Plot Figure 2c (Daylight)
egc_asv_dissim_pca_daylight = ggplot(egc_asv_dbrda_ord_points_and_meta, aes(CAP1, CAP2)) + 
  geom_point(aes(colour=Daylight, shape=Mooring_position), size=4) + 
  scale_colour_gradient2(low = "#D9D9D9", high="#611a00",
                         mid="#ff7947", midpoint=12, 
                         guide="colourbar") + 
  labs(x = "CAP1 (16%)", y = "CAP2 (7%)", colour="Daylight (h)") + 
  theme_bw() + 
  theme(legend.position="bottom", 
        legend.title = element_text(size = 14, colour="black"), 
        legend.text = element_text(size = 12, colour="black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.title.x = element_text(size = 14, colour = "black")) + 
  guides(shape=guide_legend(override.aes = list(size=5)))

### Combine three parts into single figure and export
pdf(file="Figure_2.pdf", height=5, width=14)
egc_asv_dissim_pca_aw_proportion+egc_asv_dissim_pca_ice+egc_asv_dissim_pca_daylight + 
  plot_layout(guides="collect")
dev.off()

### The legends of the plots were further repositioned in Inkscape

############################################################################################
### Figure 3 ###
############################################################################################

##########

### Figure 3. Distribution dynamics and co-occurrence of ASVs. 

##########

# Import ASV distribution dynamics and taxa info
egc_asv_dynamics_df=read.table("FRAM_RAS_EGC_ASV_dynamics_summary_and_taxa_complete.txt", header = TRUE, 
                               sep = "\t",as.is=TRUE)

# Import ASV_relative_abundance_matrix
egc_asv_rel=read.table("FRAM_RAS_EGC_ASV_relative_filtered.txt", header = TRUE, 
                       sep = "\t",as.is=TRUE,check.names=F,row.names=1)

# Import sample metadata
egc_asv_meta=read.table("FRAM_RAS_EGC_ASV_meta_refined.txt", header = TRUE, 
                        sep = "\t",check.names=F,as.is=TRUE)


egc_asv_meta$Date <- as.Date(egc_asv_meta$Date, format="%m/%d/%Y")

Figure_3a <- ggplot(egc_asv_meta, aes(x = as.Date(Date), y = Ice_cover)) + 
  geom_area(fill="grey90") + 
  geom_line(aes(x = as.Date(Date), y = AW_proportion*100), colour = "#1338BE", size = 2) +
  geom_line(aes(x = as.Date(Date), y = Daylight*4), colour="#FFD300", size = 2) +
  geom_line(aes(x = as.Date(Date), y = ifelse(Chlorophyll_a_sensor > 0, Chlorophyll_a_sensor*150, 0)), 
            colour="#85CC6F", size = 2) + 
  scale_y_continuous(sec.axis = sec_axis(~ ./150, name = "Chlorophyll a (mg m3)")) + 
  scale_x_date(date_breaks = "3 months") + 
  labs(y = "AW proportion / Ice concentration (%)", x = "Date") + 
  facet_grid(.~Mooring_position, scales="free_x") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_blank(), 
        strip.background.x = element_rect(fill = "white", colour = "black"), 
        strip.text.x = element_text(size = 16, colour = "black"))
Figure_3a

###
# Create dataframe with relative abundance per Distribution group instead of per ASV
asv_distr_group <- egc_asv_dynamics_df %>% 
  subset(., select=c(ASV_name,Distribution_group))

distr_group_rel_abund = asv_rel %>%
  tibble::rownames_to_column(., var="ASV_name") %>%
  left_join(asv_distr_group, by="ASV_name") %>% 
  melt(., id.vars=c("ASV_name","Distribution_group"), variable.name="Sample",
       value.name="Rel_abund") %>%
  aggregate(Rel_abund~Distribution_group+Sample, data=., FUN=sum) %>% 
  rename(., RAS_id = Sample) %>% 
  left_join(egc_asv_meta, by="RAS_id")

# Reformat variables and define order ready for plotting
distr_group_rel_abund$Date <- as.Date(distr_group_rel_abund$Date, 
                                                     format="%m/%d/%Y")
distr_group_rel_abund$Mooring_position <- factor(distr_group_rel_abund$Mooring_position, 
                                                       levels=unique(c("core-EGC", "MIZ")))
distr_group_rel_abund$Distribution_group <- factor(distr_group_rel_abund$Distribution_group, 
                                                     levels=unique(c("Intermittent", "Transient", "Resident")))

### Plot Figure 3b - Dynamics of community fractions across samples

Figure_3b <- ggplot(data=distr_group_rel_abund, 
                                          aes(x=as.Date(Date), y=Rel_abund)) + 
  scale_x_date(breaks = "4 months") + 
  labs(x = "Date", y = "Relative abundance (%)", fill = "Distribution group") + 
  geom_bar(aes(fill=Distribution_group), stat="identity", position="stack", width = 8) + 
  facet_grid(.~Mooring_position, scales="free_x") + 
  scale_fill_manual(values=c("Intermittent" = "#C5C6D0",
                             "Transient" = "#232023", 
                             "Resident" = "#787276")) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=14, colour="Black"),
        axis.title.y = element_text(size=16, colour="Black"),
        axis.text.x = element_text(size=14, colour="Black", angle=45, hjust=1),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 16, colour = "black"), 
        legend.text = element_text(size = 14, colour = "black"), 
        strip.background = element_rect(fill="white"), 
        strip.text.x = element_text(size = 16, colour = "black"))
Figure_3b

pdf(file = "Figure_3.pdf", height=8, width=12)
Figure_3a/Figure_3b
dev.off()

############################################################################################
### Figure 4 ###
############################################################################################

##########

### Figure 4. Sparse partial least square regression linking community structure 
### and environmental parameters.   

##########

# Import data
egc_asv_dynamics_df=read.table("FRAM_RAS_EGC_ASV_dynamics_summary_and_taxa_complete.txt", header = TRUE, 
                               sep = "\t",as.is=TRUE, check.names = F)
egc_asv_relative=read.table("FRAM_RAS_EGC_ASV_relative_filtered.txt", header = TRUE, 
                                sep = "\t",row.names=1,as.is=TRUE, check.names = F)
egc_asv_meta=read.table("FRAM_RAS_EGC_ASV_meta_refined.txt", header = TRUE, 
                            sep = "\t",as.is=TRUE,check.names = F)

# Subset and reformat metadata
egc_asv_meta_subset = egc_asv_meta %>% 
  subset(., select=c(RAS_id,Longitude,AW_proportion,
                     Daylight,Ice_cover,Ice_cover_past,Ice_edge_distance,
                     O2_concentration)) %>% 
  arrange(., RAS_id) %>% 
  textshape::column_to_rownames(., "RAS_id")
egc_asv_meta_subset <- egc_asv_meta_subset[ order(rownames(egc_asv_meta_subset)), ]

# Filter ASV data to remove low abundant, then hellinger transform and reformat
egc_asv_relative_filt = egc_asv_relative %>% 
  OTUtable::filter_taxa(., abundance=0.1, persistence = 10) %>% 
  vegan::decostand(., method="hellinger") %>%
  t(.) %>%
  as.data.frame(.)

egc_asv_relative_filt <- egc_asv_relative_filt[ order(rownames(egc_asv_relative_filt)), ]

# Run sPLS
PLS_v2 <- spls(
  egc_asv_relative_filt,  egc_asv_meta_subset,
  ncomp = 3, 
  tol = 1e-06,
  max.iter = 100,
  near.zero.var = F,
  mode = "regression")

# Find significant ASVs
cim_v2 <- cim(
  PLS_v2,
  comp = 1:3,
  cutoff = 0.4,
  dist.method = c(
    "correlation","correlation"),
  clust.method = c(
    "average","average"),
  mapping = "XY")

# Subset significant ASVs
mat_v2 <- cim_v2$mat.cor
subset_v2 = row.names(mat_v2)
subset_v2

# Filter ASV relative abundance dataframe to only include significant ASVs
asv_spls_v2_sig <- egc_asv_relative_filt[names(egc_asv_relative_filt) %in% subset_v2] 
asv_spls_v2_sig

# Test which clustering method best represents the variation in the data
# First perform hierarchical clustering with different algorithms
asv_spls_v2_hc_complete <- hclust(Dist(
  cim_v2$mat.cor, method="correlation"),"complete")	
asv_spls_v2_hc_average <- hclust(Dist(
  cim_v2$mat.cor, method="correlation"),"average")
asv_spls_v2_hc_single <- hclust(Dist(
  cim_v2$mat.cor, method="correlation"),"single")
asv_spls_v2_hc_ward <- hclust(Dist(
  cim_v2$mat.cor, method="correlation"),"ward.D")	

# Calculate cophenetic correlation for each method
comp_coph <- cophenetic(asv_spls_v2_hc_complete)
average_coph <- cophenetic(asv_spls_v2_hc_average)
single_coph <- cophenetic(asv_spls_v2_hc_single)
ward_coph <- cophenetic(asv_spls_v2_hc_single)

# Print cophenetic values
paste("Complete cophenetic correlation: ", cor(Dist(
  cim_v2$mat.cor), comp_coph))
paste("Average cophenetic correlation: ", cor(Dist(
  cim_v2$mat.cor), average_coph))
paste("Single cophenetic correlation: ", cor(Dist(
  cim_v2$mat.cor), single_coph))
paste("Ward cophenetic correlation: ", cor(Dist(
  cim_v2$mat.cor), ward_coph))

### The highest cophenetic correlation was observed using the 'average' method

# Determine the optimal number of clusters that best represents variation
# between the ASVs

asv_spls_v2_hc_elbow <- 
  fviz_nbclust(asv_spls_v2_sig,
               diss = Dist(
                 cim_v2$mat.cor, method="correlation"), hcut, method = "wss") + 
  labs(title = "Hierarchical clust UPGMA - Elbow method")
asv_spls_v2_hc_elbow

asv_spls_v2_hc_gap_statistics <- 
  fviz_nbclust(asv_spls_v2_sig, diss = Dist(
    cim_v2$mat.cor, method="correlation"), hcut,
               nstart = 25, method = "gap_stat",
               nboot = 100) + 
  labs(subtitle = "Gap statistic method")
asv_spls_v2_hc_gap_statistics

### Based on these two calculations, we could determine that between 7 - 10
### clusters is likely optimal

# Repeat sPLS analysis using only the significant ASVs along with the 
# average clustering method

# Final PLS with subset
PLS_v2_final <- spls(
  asv_spls_v2_sig, egc_asv_meta_subset, 
  ncomp = 3, 
  tol = 1e-06,
  max.iter = 100,
  near.zero.var = F,
  mode = "regression")

cim_v2_final <- cim(
  PLS_v2_final,
  comp = 1:3,
  dist.method = c(
    "correlation","correlation"), 
  clust.method=c(
    "average", "average"), 
  mapping = "XY",
  color = colorRampPalette(c(
    "#543005","#ffffff","#003C30"))(40))

# Extract correlations
asv_spls_sig_v2_cor_mat = data.frame(cim_v2_final$mat.cor) %>%
  mutate(across(everything(), round, 3))
dim(asv_spls_sig_v2_cor_mat)

# Compute hierarchical clustering
asv_spls_v2_hc <- hclust(Dist(
  cim_v2_final$mat.cor, method="correlation"),"average")	

# Visualise
plot(asv_spls_v2_hc)

# Define 7 clusters and plot on dendrogram
asv_spls_v2_cutree_7 <- cutree(asv_spls_v2_hc, k=7)
plot(as.phylo(asv_spls_v2_hc), tip.color=asv_spls_v2_cutree_7, cex=0.2)

# Define 8 clusters and plot on dendrogram
asv_spls_v2_cutree_8 <- cutree(asv_spls_v2_hc, k=8)
plot(as.phylo(asv_spls_v2_hc), tip.color=asv_spls_v2_cutree_8, cex=0.2)

# Define 9 clusters and plot on dendrogram
asv_spls_v2_cutree_9 <- cutree(asv_spls_v2_hc, k=9)
plot(as.phylo(asv_spls_v2_hc), tip.color=asv_spls_v2_cutree_9, cex=0.2)

### Using 9 clusters appears to best separate out the distinct branches on
### the dendrogram.

# Obtain taxonomic affiliation of
# Extract the asv assignment in the 9 clusters and combine with other
# ASV data
clusters <- cbind(asv_spls_sig_v2_cor_mat, asv_spls_v2_cutree_9) %>%
  tibble::rownames_to_column("asv") %>%
  mutate(cluster=paste0("C", asv_spls_v2_cutree_9)) %>% 
  dplyr::slice(mixedorder(asv))

# Define cluster colors
colors <- c(
  "C1" = "#db6d00",
  "C2" = "#920000",
  "C3" = "#b6dbff",
  "C4" = "#006ddb",
  "C5" = "#009292",
  "C6" = "#ffb6db",
  "C7" = "#ffff6d",
  "C8" = "#000000",
  "C9" = "#004949")

# Add cluster colors 
clusters$color <- colors[match(
  clusters$cluster, names(colors))]

# Plot the result with the defined cluster annotations  
cim(
  PLS_v2_final,
  dist.method = c(
    "correlation","correlation"), 
  clust.method=c(
    "average","average"), 
  margins = c(7,20), 
  row.cex = 0.4,
  row.names = F,
  mapping = "XY",
  row.sideColors = clusters$color,
  legend = list(
    legend = names(colors),
    col = c(colors),
    cex = 1.44,
    title = "Cluster"), 
  color = colorRampPalette(c(
    "#000000", "#ffffff", "#BDA55D"))(40))
  #save = "pdf",
  #name.save = "Figure_4a_new") 

### Inspecting the result - 
### Clusters C3 and C4 only include a few ASVs and largely both show
### strong negative correlations to daylight with only slight variations in 
### correlation to water mass.
### As a result, to reduce the number of small clusters and increase the focus
### on distinct conditions, these clusters will be grouped together 

clusters_and_correlation_mat_final <- cbind(asv_spls_sig_v2_cor_mat, asv_spls_v2_cutree_9) %>%
  tibble::rownames_to_column("asv") %>%
  mutate(cluster=paste0("C", asv_spls_v2_cutree_9)) %>% 
  dplyr::slice(mixedorder(asv)) %>% 
  mutate(cluster=case_when(
    cluster %in% c("C3","C4")~"C3",
    cluster =="C5"~"C4",
    cluster =="C6"~"C5",
    cluster =="C7"~"C6",
    cluster =="C8"~"C7",
    cluster=="C9"~"C8",
    TRUE~cluster))

# Define cluster colors
colors <- c(
  "C1" = "#db6d00",
  "C2" = "#920000",
  "C3" = "#b6dbff",
  "C4" = "#006ddb",
  "C5" = "#009292",
  "C6" = "#ffb6db",
  "C7" = "#ffff6d",
  "C8" = "#000000")

# Add cluster colors 
clusters_and_correlation_mat_final$color <- colors[match(
  clusters_and_correlation_mat_final$cluster, names(colors))]

# Plot the final sPLS matrix result with the defined cluster annotations  
cim(
  PLS_v2_final,
  dist.method = c(
    "correlation","correlation"), 
  clust.method=c(
    "average","average"), 
  margins = c(7,20), 
  row.cex = 0.4,
  row.names = F,
  mapping = "XY",
  row.sideColors = clusters_and_correlation_mat_final$color,
  legend = list(
    legend = names(colors),
    col = c(colors),
    cex = 1.44,
    title = "Cluster"), 
  color = colorRampPalette(c(
    "#000000", "#ffffff", "#BDA55D"))(40),
  save = "pdf",
  name.save = "Figure_4a_new") 

# Combine cluster information with previous asv dynamics summary dataframe
# and export
egc_asv_cluster_temp = clusters_and_correlation_mat_final %>%
  subset(., select=c(asv,cluster)) %>%
  rename(., ASV_name = asv)

egc_asv_dynamics_df_with_cluster <- left_join(egc_asv_dynamics_df, egc_asv_cluster_temp, by="ASV_name")

write.table(egc_asv_dynamics_df_with_cluster, 
            file="FRAM_RAS_EGC_ASV_dynamics_summary_and_taxa_complete_with_spls_cluster.txt",
            sep="\t")

#####

### Figure 4b

# Create dataframe with only most prominent ASVs in each cluster 
# Then group into genera

egc_asv_dynamics_df_with_cluster=read.table("FRAM_RAS_EGC_ASV_dynamics_summary_and_taxa_complete_with_spls_cluster.txt", header = TRUE, 
                               sep = "\t",as.is=TRUE,check.names=F)
egc_asv_relative=read.table("FRAM_RAS_EGC_ASV_relative_filtered.txt", header = TRUE, 
                            sep = "\t",as.is=TRUE,check.names=F)

egc_asv_relative_mod = egc_asv_relative %>% 
  tibble::rownames_to_column(., var="ASV_name")

egc_asv_most_prominent_in_each_cluster_summary_df= egc_asv_dynamics_df_with_cluster %>% 
  filter(cluster != "NA") %>% 
  filter(Max_rel_abund > 1) %>% 
  subset(., select=c(ASV_name,Class,Genus,cluster)) %>% 
  left_join(egc_asv_relative_mod, by="ASV_name") %>% 
  reshape2::melt(., id.vars=c("ASV_name", "Class", "Genus", "cluster"), variable.name="Sample",
       value.name="rel_abund") %>% 
  aggregate(rel_abund~cluster+Class+Genus+Sample, data=., FUN=sum) %>% 
  arrange(., Genus,desc(rel_abund)) %>% 
  group_by(cluster,Genus) %>% 
  slice(1) %>% 
  arrange(., cluster,Genus)


# Plot Figure 4b

pls_cluster_colours <- c(
  "C1" = "#db6d00",
  "C2" = "#920000",
  "C3" = "#b6dbff",
  "C4" = "#006ddb",
  "C5" = "#009292",
  "C6" = "#ffb6db",
  "C7" = "#ffff6d",
  "C8" = "#000000")

egc_asv_most_prominent_in_each_cluster_summary_df$cluster <- 
  factor(egc_asv_most_prominent_in_each_cluster_summary_df$cluster,
         levels = unique(egc_asv_most_prominent_in_each_cluster_summary_df$cluster))
egc_asv_most_prominent_in_each_cluster_summary_df$Genus <- 
  factor(egc_asv_most_prominent_in_each_cluster_summary_df$Genus,
         levels = unique(egc_asv_most_prominent_in_each_cluster_summary_df$Genus))

pdf(file="Figure_4b_new.pdf", height=10, width=8)
ggplot(data=egc_asv_most_prominent_in_each_cluster_summary_df) + 
  geom_point(aes(x=cluster, y=Genus, size=rel_abund, colour=cluster), stat="identity",
             position="identity") + 
  scale_colour_manual(values=pls_cluster_colours) + 
  scale_size_continuous(range = c(2, 8), breaks=c(2,4,8,16), 
                        name="Maximum relative abundance (%)") + 
  scale_y_discrete(limits=rev) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black")) + 
  guides(colour=guide_legend(override.aes = list(size=5)))
dev.off()

#####

### Figures 4a and 4b were both exported and combined in Inkscape

############################################################################################
### Figure 5 ###
############################################################################################

##########

### Figure 5. Ribosomal protein-based phylogenetic tree of MAGs from this study 
### and those previously recovered from the Fram Strait. 

##########

### The figure was created in iToL


############################################################################################
### Figure 6 ###
############################################################################################

##########

### Figure 6. Temporal dynamics of signature populations. 

##########

# Import data
egc_asv_meta=read.table("FRAM_RAS_EGC_ASV_meta_refined.txt", header = TRUE, 
                        sep = "\t",as.is=TRUE,check.names=F)
egc_asv_meta$Date <- as.Date(egc_asv_meta$Date, format="%m/%d/%Y")

# Plot key metadata variables over temporal scales

egc_asv_meta_summary_plot <- ggplot(egc_asv_meta, aes(x = as.Date(Date), y = Ice_cover)) + 
  geom_area(fill="grey90") + 
  geom_line(aes(x = as.Date(Date), y = AW_proportion*100), colour = "#1338BE", size = 2) +
  geom_line(aes(x = as.Date(Date), y = Daylight*4), colour="#FFD300", size = 2) +
  geom_line(aes(x = as.Date(Date), y = ifelse(Chlorophyll_a_sensor > 0, Chlorophyll_a_sensor*150, 0)), colour="#85CC6F", size = 2) + 
  scale_y_continuous(sec.axis = sec_axis(~ ./150, name = "Chlorophyll a (mg m3)")) + 
  scale_x_date(date_breaks = "3 months") + 
  labs(y = "AW proportion / Ice concentration (%)", x = "Date") + 
  facet_grid(.~Mooring_position, scales="free_x") + 
  theme_bw() + 
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_blank(), 
        strip.background.x = element_rect(fill = "white", colour = "black"), 
        strip.text.x = element_text(size = 16, colour = "black"))
egc_asv_meta_summary_plot

### 

# Import ASV data

egc_asv_dynamics_df_with_cluster=read.table("FRAM_RAS_EGC_ASV_dynamics_summary_and_taxa_complete_with_spls_cluster.txt", header = TRUE, 
                                            sep = "\t",as.is=TRUE,check.names=F)
egc_asv_relative=read.table("FRAM_RAS_EGC_ASV_relative_filtered.txt", header = TRUE, 
                            sep = "\t",as.is=TRUE,check.names=F)

# Import ASV to MAG information

egc_asv_to_mag=read.table("FRAM_RAS_EGC_ASV_to_MAG.txt", header = TRUE, 
                            sep = "\t",as.is=TRUE,check.names=F)

egc_asv_dynamics_df_with_cluster_and_signatures = egc_asv_dynamics_df_with_cluster %>% 
  left_join(egc_asv_to_mag, by="ASV_name") %>% 
  mutate(Signature_population = case_when(
    MAG_name != "NA" & cluster != "NA" ~ "Signature"
  ))

# Export this dataframe as it now contains all information used throughout the 
# analysis

write.table(egc_asv_dynamics_df_with_cluster_and_signatures,
             file="FRAM_RAS_EGC_ASV_complete_summary_information.txt",
             sep="\t")

# Reformat ASV relative abundance matrix

egc_asv_relative_mod = egc_asv_relative %>% 
  tibble::rownames_to_column(., var="ASV_name")

# Subset information from metadata table that will be used

egc_asv_meta_mod <- egc_asv_meta %>% 
  subset(., select=c(RAS_id,Mooring_position,Date)) %>%
  rename(., Sample = RAS_id)

# Filter ASV information to only include signature populations and combine
# with metadata

egc_asv_pls_signatures_with_meta = egc_asv_dynamics_df_with_cluster_and_signatures %>%
  filter(Signature_population == "Signature") %>% 
  subset(., select=c(ASV_name, cluster, Genus)) %>% 
  mutate(Name_and_genus=paste(ASV_name, Genus, sep="-")) %>% 
  left_join(egc_asv_relative_mod, by="ASV_name") %>% 
  reshape2::melt(id.vars=c("ASV_name","cluster","Genus","Name_and_genus"), variable.name="Sample",
       value.name="rel_abund") %>%
  left_join(egc_asv_meta_mod, by="Sample")
View(egc_asv_pls_signatures_with_meta)

# Define factor levels

egc_asv_pls_signatures_with_meta$Mooring_position <- factor(egc_asv_pls_signatures_with_meta$Mooring_position, 
                                                   levels=unique(c("core-EGC","MIZ")))
egc_asv_pls_signatures_with_meta$Date <- as.Date(egc_asv_pls_signatures_with_meta$Date, format="%m/%d/%Y")

### Plot
# Plot sPLS1 signature populations temporal dynamics

pls_c1_temporal= egc_asv_pls_signatures_with_meta %>% 
  filter(cluster == "C1") %>% 
  ggplot(data=., aes(x=as.Date(Date), y=rel_abund)) + 
  scale_x_date(breaks = "4 months") + 
  labs(x = "Date", y = "Relative abundance (%)", 
       fill = "C1 signature populations") + 
  geom_bar(aes(fill=Name_and_genus), stat="identity", position="stack", width = 8) + 
  scale_fill_manual(values = c("#A6CEE3", "#1F78B4")) + 
  facet_grid(.~Mooring_position, scales="free_x") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12, colour="Black"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 10, colour = "black"), 
        legend.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        plot.background = element_blank()) + 
  guides(fill=guide_legend(ncol=2))

# Plot sPLS2 signature populations temporal dynamics

pls_c2_temporal = egc_asv_pls_signatures_with_meta %>% 
  filter(cluster == "C2") %>% 
  ggplot(data=., aes(x=as.Date(Date), y=rel_abund)) + 
  scale_x_date(breaks = "4 months") + 
  labs(x = "Date", y = "Relative abundance (%)", 
       fill = "C2 signature populations") + 
  geom_bar(aes(fill=Name_and_genus), stat="identity", position="stack", width = 8) + 
  scale_fill_manual(values = "#B2DF8A") + 
  facet_grid(.~Mooring_position, scales="free_x") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12, colour="Black"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 10, colour = "black"), 
        legend.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        plot.background = element_blank())

# Plot sPLS3 signature populations temporal dynamics

pls_c3_temporal = egc_asv_pls_signatures_with_meta %>% 
  filter(cluster == "C3") %>% 
  ggplot(data=., aes(x=as.Date(Date), y=rel_abund)) + 
  scale_x_date(breaks = "4 months") + 
  labs(x = "Date", y = "Relative abundance (%)", 
       fill = "C3 signature populations") + 
  geom_bar(aes(fill=Name_and_genus), stat="identity", position="stack", width = 8) + 
  scale_fill_manual(values = c("#33A02C","#FB9A99")) + 
  facet_grid(.~Mooring_position, scales="free_x") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12, colour="Black"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 10, colour = "black"), 
        legend.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        plot.background = element_blank())+ 
  guides(fill=guide_legend(ncol=2))

# Plot sPLS4 signature populations temporal dynamics

pls_c4_temporal = egc_asv_pls_signatures_with_meta %>% 
  filter(cluster == "C4") %>% 
  ggplot(data=., aes(x=as.Date(Date), y=rel_abund)) + 
  scale_x_date(breaks = "4 months") + 
  labs(x = "Date", y = "Relative abundance (%)", 
       fill = "C4 signature populations") + 
  geom_bar(aes(fill=Name_and_genus), stat="identity", position="stack", width = 8) + 
  scale_fill_manual(values = c("#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A",
                    "#B15928")) + 
  facet_grid(.~Mooring_position, scales="free_x") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12, colour="Black"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 10, colour = "black"), 
        legend.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        plot.background = element_blank())+ 
  guides(fill=guide_legend(ncol=3))

# Plot sPLS6 signature populations temporal dynamics

pls_c6_temporal = egc_asv_pls_signatures_with_meta %>% 
  filter(cluster == "C6") %>% 
  ggplot(data=., aes(x=as.Date(Date), y=rel_abund)) + 
  scale_x_date(breaks = "4 months") + 
  labs(x = "Date", y = "Relative abundance (%)", 
       fill = "C6 signature populations") + 
  geom_bar(aes(fill=Name_and_genus), stat="identity", position="stack", width = 8) + 
  scale_fill_manual(values = "#A6CEE3") +  
  facet_grid(.~Mooring_position, scales="free_x") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12, colour="Black"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 10, colour = "black"), 
        legend.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        plot.background = element_blank())

# Plot sPLS7 signature populations temporal dynamics

pls_c7_temporal = egc_asv_pls_signatures_with_meta %>% 
  filter(cluster == "C7") %>% 
  ggplot(data=., aes(x=as.Date(Date), y=rel_abund)) + 
  scale_x_date(breaks = "4 months") + 
  labs(x = "Date", y = "Relative abundance (%)", 
       fill = "C7 signature populations") + 
  geom_bar(aes(fill=Name_and_genus), stat="identity", position="stack", width = 8) + 
  scale_fill_manual(values = "#1F78B4") +  
  facet_grid(.~Mooring_position, scales="free_x") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12, colour="Black"),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 10, colour = "black"), 
        legend.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        plot.background = element_blank())

# Plot sPLS8 signature populations temporal dynamics

pls_c8_temporal = egc_asv_pls_signatures_with_meta %>% 
  filter(cluster == "C8") %>% 
  ggplot(data=., aes(x=as.Date(Date), y=rel_abund)) + 
  scale_x_date(breaks = "4 months") + 
  labs(x = "Date", y = "Relative abundance (%)", 
       fill = "C8 signature populations") + 
  geom_bar(aes(fill=Name_and_genus), stat="identity", position="stack", width = 8) + 
  scale_fill_manual(values = c("#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F")) +
  facet_grid(.~Mooring_position, scales="free_x") + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12, colour="Black"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, colour="black", angle=45, hjust=1),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 12, colour = "black"), 
        legend.text = element_text(size = 10, colour = "black"), 
        legend.background = element_blank(),
        strip.background = element_blank(), 
        strip.text.x = element_blank(),
        plot.background = element_blank()) +
  guides(fill=guide_legend(ncol=3))

### Plot Figure 6
### Combine all plots into single figure

library(patchwork)
pdf(file="Figure_6_new.pdf", height=12, width=14)
(egc_asv_meta_summary_plot/pls_c1_temporal/pls_c2_temporal/
  pls_c3_temporal/pls_c4_temporal/pls_c6_temporal/pls_c7_temporal/pls_c8_temporal)+
  plot_layout(guides="collect")
dev.off()

### Plot was further modified in Inkscape

############################################################################################
### Figure 7 ###
############################################################################################

##########

### Figure 7. Clustering of sample functional gene composition based on Bray-Curtis dissimilarities.

##########

# Import pre-normalised functional profile matrix

egc_read_func_long=read.table("FRAM_RAS_EGC_community_functional_profile_norm_long.txt", 
                              header = TRUE, sep = "\t",as.is=TRUE,check.names=F)
egc_meta=read.table("FRAM_RAS_EGC_MG_meta.txt", 
                           header = TRUE, sep = "\t",as.is=TRUE,check.names=F)

# Reformat dataframe

egc_read_func_wide = egc_read_func_long %>% 
  aggregate(Norm_count~Gene+Sample, data=., FUN = sum) %>% 
  reshape2::dcast(Gene~Sample, value.var = "Norm_count", data=.) %>% 
  tibble::column_to_rownames(., var="Gene")

egc_read_func_wide[is.na(egc_read_func_wide)]<-0

# Filter out very low abundant functions

egc_read_func_wide_filt <- OTUtable::filter_taxa(egc_read_func_wide, abundance=0.001, persistence = 20)
dim(egc_read_func_wide_filt)

# Compute dissimilarity matrix

egc_read_func_wide_filt_t <- t(egc_read_func_wide_filt)
egc_read_func_wide_bray <- vegdist(egc_read_func_wide_filt_t, method="bray", upper=FALSE)

# Perform hierarchical clustering using different algorithms

egc_read_func_wide_clust_avg <- hclust(egc_read_func_wide_bray, method="average")
egc_read_func_wide_clust_comp <- hclust(egc_read_func_wide_bray, method="complete")
egc_read_func_wide_clust_ward <- hclust(egc_read_func_wide_bray, method="ward.D")

par(mfrow=c(3,1))
plot(egc_read_func_wide_clust_avg)
plot(egc_read_func_wide_clust_comp)
plot(egc_read_func_wide_clust_ward)

# Although slight variation observed, the same two distinct clusters were
# consistently formed

# Check which method is likely better at capturing dissimilarity

# Calculate cophenetic correlation for each method
egc_func_avg_coph <- cophenetic(egc_read_func_wide_clust_avg)
egc_func_comp_coph <- cophenetic(egc_read_func_wide_clust_comp)
egc_func_ward_coph <- cophenetic(egc_read_func_wide_clust_ward)

# Print cophenetic values
paste("Complete cophenetic correlation: ", cor(egc_read_func_wide_bray,
                                               egc_func_comp_coph))
paste("Average cophenetic correlation: ", cor(egc_read_func_wide_bray, 
                                              egc_func_avg_coph))
paste("Ward cophenetic correlation: ", cor(egc_read_func_wide_bray, 
                                           egc_func_ward_coph))

# Based on the highest cophenetic correlation, the average algorithm will be used

# Define clusters and combine with metadata

egc_func_clusters <- cutree(egc_read_func_wide_clust_avg, k=2)

# Create dataframe with cluster information
egc_func_clusters_df = egc_func_clusters %>% 
  as.data.frame() %>%
  rename(., Clusters = ".") %>% 
  tibble::rownames_to_column(., "RAS_id")

# Save cluster assignment output for use in Figure 8/9/10
write.table(egc_func_clusters_df, file="FRAM_RAS_EGC_MG_cluster_assignment.txt", 
            sep="\t")

# Combine cluster information with the main metadata variables
egc_func_metadata = egc_func_clusters_df %>%
  left_join(egc_meta, by="RAS_id") %>% 
  subset(., select=c(RAS_id,Ice_cover,Ice_edge_distance,
                     Ice_cover_past,Daylight,AW_proportion,O2_concentration)) %>%
  tibble::column_to_rownames(., "RAS_id")

# Check for colinearity between metadata variables
ggpairs(egc_func_metadata)

### Ice edge distance and ice cover highly colinear
### Ice cover and O2 concentration highly colinear
### Ice cover will be kept and the others removed

# Standardize metadata
egc_func_meta_stand_with_clusters = egc_func_metadata %>% 
  subset(., select=-c(O2_concentration,Ice_edge_distance)) %>%
  decostand(., method="standardize", MARGIN=2) %>%
  tibble::rownames_to_column(., "RAS_id") %>%
  left_join(egc_func_clusters_df, by="RAS_id")

# Run MANOVA to determine which metadata variables are significant between
# clusters 

depend_vars <- cbind(egc_func_meta_stand_with_clusters$Ice_cover,
                     egc_func_meta_stand_with_clusters$Ice_cover_past,
                     egc_func_meta_stand_with_clusters$Daylight,
                     egc_func_meta_stand_with_clusters$AW_proportion)

egc_read_func_meta_cluster_manova <- manova(depend_vars~Clusters, data=egc_func_meta_stand_with_clusters)
summary(egc_read_func_meta_cluster_manova)
summary.aov(egc_read_func_meta_cluster_manova)

### This identified Ice_cover as the only significant variable

### Prepare metadata dataframe for plotting
egc_func_metadata_for_plotting = egc_func_metadata %>% 
  subset(., select=-c(O2_concentration,Ice_edge_distance)) %>%
  tibble::rownames_to_column(., "RAS_id") %>%
  left_join(egc_func_clusters_df, by="RAS_id") %>%
  mutate(AW_proportion_2 = AW_proportion*100) %>% 
  subset(., select=-c(AW_proportion)) %>% 
  rename(., AW_proportion = AW_proportion_2)

# Sort by sample ID as this has to match the dissimilarity matrix order
egc_func_metadata_for_plotting <- egc_func_metadata_for_plotting[order(row.names(egc_func_metadata_for_plotting)),]

# Define clusters as a clustering object
egc_read_func_wide_bray <- vegdist(egc_read_func_wide_filt_t, method="bray", upper=FALSE)
egc_read_func_clusters <- hcut(egc_read_func_wide_bray, k=2, hc_method="average", isdiss=TRUE)

# Define dissimilarity matrix as a matrix
egc_read_func_wide_bray <- as.matrix(egc_read_func_wide_bray)
egc_read_func_wide_bray

### Plot dissimilarities at heatmap with metadata information also
taxa_heatmap_colours <- colorRamp2(c(0.05, 0.3), c("#150c25", "#FFFFE5"))

### Barplot annotation
column_ha_bars = HeatmapAnnotation("Ice coverage (%)" = anno_barplot(egc_func_metadata_for_plotting$Ice_cover, 
                                                                     axis_param = list(at = c(25,50,75,100))), 
                                   "Atlantic water proportion (%)" = anno_barplot(egc_func_metadata_for_plotting$AW_proportion, 
                                                                                  axis_param = list(at = c(25,50,75,100))), 
                                   "Daylight (h)" = anno_barplot(egc_func_metadata_for_plotting$Daylight, 
                                                                 axis_param = list(at = c(6,12,18,24))),
                                   height = unit(6, "cm"),
                                   annotation_name_side = "left",
                                   annotation_name_rot = 1)

#
egc_read_func_heatmap <- Heatmap(egc_read_func_wide_bray, 
                                 #rect_gp = gpar(type = "none"), 
                                 column_dend_side = "bottom",
                                 col = taxa_heatmap_colours,
                                 show_column_names = FALSE,
                                 top_annotation=column_ha_bars,
                                 #bottom_annotation=col_cluster_symbol,
                                 #left_annotation=row_cluster_symbol,
                                 row_names_side = c("left"), 
                                 cluster_columns=egc_read_func_clusters,
                                 cluster_rows=egc_read_func_clusters,
                                 clustering_method_columns = "average",
                                 clustering_method_rows = "average", 
                                 column_dend_height = unit(60, "pt"),
                                 heatmap_legend_param = list(title="Bray Curtis dissimilarity", 
                                                             legend_height = unit(4, "cm"), 
                                                             labels_gp=gpar(fontsize=12), 
                                                             title_gp=gpar(fotnsize=14)))

egc_read_func_heatmap

pdf(file="Figure_7.pdf", height=7, width=9)
egc_read_func_heatmap
dev.off()


############################################################################################
### Figure 8 and Figure 9 ###
############################################################################################

### Both figure 8 and figure 9 show different resolutions of significantly
### enriched genes derived from the same analysis
### As such, the first section of the next code will perform the analysis
### and the respective figures are generated at the end

#####

### Aldex2 analysis

#####

# Import data
egc_read_func_long=read.table("FRAM_RAS_EGC_community_functional_profile_norm_long.txt", 
                              header = TRUE, sep = "\t",as.is=TRUE,check.names=F)
egc_meta=read.table("FRAM_RAS_EGC_MG_meta.txt", 
                    header = TRUE, sep = "\t",as.is=TRUE,check.names=F)
egc_mg_cluster_assignment=read.table("FRAM_RAS_EGC_MG_cluster_assignment.txt", 
                    header = TRUE, sep = "\t",as.is=TRUE,check.names=F)

# Reformat

egc_read_func_wide = egc_read_func_long %>% 
  aggregate(Norm_count~Gene+Sample, data=., FUN = sum) %>% 
  reshape2::dcast(Gene~Sample, value.var = "Norm_count", data=.) %>% 
  tibble::column_to_rownames(., var="Gene")

egc_read_func_wide[is.na(egc_read_func_wide)]<-0

# Filter out very low abundant functions and reformat dataframe to long format

egc_read_func_long_filt <- egc_read_func_wide %>% 
  OTUtable::filter_taxa(., abundance=0.001, persistence = 20) %>% 
  tibble::rownames_to_column(., "Gene") %>% 
  reshape2::melt(., id.vars=c("Gene"), 
       variable.name="RAS_id", value.name="Norm_count")
View(egc_read_func_long_filt)

# Combine with metadata and perform additional reformating of variable names
egc_mg_clusters_tidy = egc_mg_cluster_assignment %>% 
  mutate(Cluster = case_when(
    Clusters == 1 ~ "High-ice cluster",
    Clusters == 2 ~ "Low-ice cluster"
  )) %>% 
  subset(., select=-c(Clusters)) %>% 
  arrange(., RAS_id)

egc_read_func_long_filt_with_cluster = egc_read_func_long_filt %>% 
  left_join(egc_mg_cluster_assignment, by="RAS_id") %>% 
  mutate(Cluster = case_when(
    Clusters == 1 ~ "High-ice cluster",
    Clusters == 2 ~ "Low-ice cluster"
  )) %>% 
  subset(., select=-c(Clusters)) %>% 
  arrange(., RAS_id)

functional_matrix_for_aldex <- egc_read_func_long_filt %>% 
  mutate_at(c("Norm_count"), as.numeric) %>%
  mutate(count=Norm_count*1000) %>% #scaling up values as absolute integers needed)
  mutate_at(c("count"), round, 0) %>% 
  reshape2::dcast(Gene~RAS_id, value.var="count", data=.,) %>% 
  column_to_rownames("Gene")

View(functional_matrix_for_aldex)
View(egc_mg_clusters_tidy)

# Run Aldex2
egc_mg_func_aldex <- aldex(
  functional_matrix_for_aldex,
  factor(egc_mg_clusters_tidy$Cluster), 
  denom = "zero",
  test = "kw",
  mc.samples = 128)
egc_mg_func_aldex

# select significant results
egc_mg_func_aldex_sig_genes <- rownames(egc_mg_func_aldex)[
  egc_mg_func_aldex$glm.eBH < 0.05 & egc_mg_func_aldex$kw.ep < 0.05]
View(egc_mg_func_aldex_sig_genes)

# Create matrix for significant hits, combine with metadata and calculate CLR
egc_mg_func_aldex_sig_genes_clr = subset(
  functional_matrix_for_aldex, rownames(functional_matrix_for_aldex) %in% egc_mg_func_aldex_sig_genes) %>% 
  rgr:: clr() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., var="Gene") %>%
  reshape2::melt(id.vars="Gene", variable.name="RAS_id", value.name="CLR_values") %>%
  left_join(egc_mg_clusters_tidy, by="RAS_id")
View(egc_mg_func_aldex_sig_genes_clr)

# Determine which group genes are enriched under and add to dataframe
egc_mg_func_gene_enrichment_category = egc_mg_func_aldex_sig_genes_clr %>% 
  aggregate(CLR_values~Gene+Cluster, data=., FUN=mean) %>% 
  reshape2::dcast(., Gene~Cluster, value.var="CLR_values") %>%
  mutate(Enrichment = case_when(
    `High-ice cluster` > `Low-ice cluster` ~ "HI_enriched",
    `Low-ice cluster` > `High-ice cluster` ~ "LI_enriched"
  )) %>% 
  subset(., select=c(Gene,Enrichment))

egc_mg_func_aldex_sig_genes_clr_and_group_enrichment = 
  left_join(egc_mg_func_aldex_sig_genes_clr,
            egc_mg_func_gene_enrichment_category,
            by="Gene")

# Import gene names for KEGG annotations
kegg_info = read.table("KEGG_ID_to_description.csv", header = TRUE, 
                       sep = ",",as.is=TRUE, check.names=F)

# Add KO descriptions to dataframe
egc_mg_func_aldex_sig_genes_clr_and_group_enrichment_updated = egc_mg_func_aldex_sig_genes_clr_and_group_enrichment %>% 
  mutate(description = "NA") %>% 
  left_join(kegg_info, by="Gene") %>% 
  subset(., select=-c(description.x)) %>% 
  mutate(KEGG_description = description.y) %>% 
  subset(., select=-c(description.y))
View(egc_mg_func_aldex_sig_genes_clr_and_group_enrichment_updated)

# Export dataframe
write.table(egc_mg_func_aldex_sig_genes_clr_and_group_enrichment_updated, 
            file="FRAM_RAS_EGC_community_functional_profile_ALDEX_sig_enr_genes_long.txt", 
            sep = "\t")

#####

### At this point, the output table produced above was manually inspected,
### Gene-encoding functions were checked and two subset dataframes were produced 
### that were used to make the following figures

##########

### Figure 8. Summary of substrate uptake and degradation-related genes 
### enriched with high- and low-ice coverage. .

##########

# Import grouped significant gene hits
aldex_enriched_genes_grouped =read.table("FRAM_RAS_EGC_community_functional_profile_ALDEX_sig_enr_genes_grouped.csv", header = TRUE, 
                                                    sep = ",",as.is=TRUE)

aldex_enriched_genes_grouped$Metabolism <- factor(aldex_enriched_genes_grouped$Metabolism, levels = unique(aldex_enriched_genes_grouped$Metabolism))
aldex_enriched_genes_grouped$Metabolism_group <- factor(aldex_enriched_genes_grouped$Metabolism_group, levels = unique(aldex_enriched_genes_grouped$Metabolism_group))

# Create boxplot
aldex_enriched_genes_grouped_boxplot <- 
  aldex_enriched_genes_grouped %>% 
  filter(Metabolism_group == "Substrates") %>% 
  ggplot(.,) +
  geom_boxplot(aes(
    x=Metabolism, y=CLR_values, 
    fill=cluster), outlier.shape = NA) + 
  geom_hline(
    yintercept=0, size=0.5,
    linetype="dashed", 
    color="gray33") + 
  scale_fill_manual(values=c(
    "High-ice cluster"="darkcyan",
    "Low-ice cluster"="lightgoldenrod2")) +  
  labs(y = "CLR-transformed normalised gene count", fill = "Ice coverage") +
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        legend.position = "right",
        axis.title.y = element_text(size = 14, colour = "black", vjust=20), 
        axis.text.x = element_text(size = 11, colour = "black", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 14, colour = "black"),
        panel.grid.minor=element_blank(), 
        plot.margin = margin(0.1,0.1,0.1,2.5,"cm"))
aldex_enriched_genes_grouped_boxplot

pdf(file="Figure_8.pdf", height = 8, width = 12)
aldex_enriched_genes_grouped_boxplot
dev.off()

##########

### Figure 9. Figure 9. Selected genes involved in the uptake and degradation 
### of organic and inorganic compounds enriched under high- and low-ice conditions. 

##########

# Import data
aldex_enriched_genes_detailed =read.table("FRAM_RAS_EGC_community_functional_profile_ALDEX_sig_enr_genes_detailed.csv", header = TRUE, 
                                                        sep = ",",as.is=TRUE)

aldex_enriched_genes_detailed$Metabolism <- factor(aldex_enriched_genes_detailed$Metabolism, levels = unique(aldex_enriched_genes_detailed$Metabolism))
aldex_enriched_genes_detailed$Gene_short <- factor(aldex_enriched_genes_detailed$Gene_short, levels = unique(aldex_enriched_genes_detailed$Gene_short))

aldex_enriched_genes_expanded_boxplot <- 
  ggplot(data=aldex_enriched_genes_detailed) +
  geom_boxplot(aes(
    x=Gene_short, y=CLR_value, 
    fill=cluster), outlier.shape = NA) + 
  geom_hline(
    yintercept=0, size=0.5,
    linetype="dashed", 
    color="gray33") + 
  scale_fill_manual(values=c(
    "High-ice cluster"="darkcyan",
    "Low-ice cluster"="lightgoldenrod2")) +  
  labs(y = "CLR-transformed normalised gene count", fill = "Ice coverage") +
  facet_wrap(~Metabolism, scales="free")+
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        legend.position = "right",
        axis.title.y = element_text(size = 14, colour = "black", vjust=20), 
        axis.text.x = element_text(size = 11, colour = "black", angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 14, colour = "black"),
        panel.grid.minor=element_blank(), 
        plot.margin = margin(0.1,0.1,0.1,2.5,"cm"))
aldex_enriched_genes_expanded_boxplot

pdf(file="Figure_9.pdf", height = 12, width = 12)
aldex_enriched_genes_expanded_boxplot
dev.off()



###################

###################

### Supplementary Figure plotting

# Import ASV distribution dynamics and taxa info
egc_asv_dynamics_df=read.table("FRAM_RAS_EGC_ASV_dynamics_summary_and_taxa_complete.txt", header = TRUE, 
                               sep = "\t",as.is=TRUE)
###
# Create dataframe with distrbution group and ASV distribution dynamics 

asv_max_abund_moorings = egc_asv_dynamics_df %>%
  subset(., select=c(ASV_name,Mooring_presence,Max_rel_abund,Num_samples_present,
                     Distribution_group))

asv_max_abund_moorings$Mooring_presence <- factor(asv_max_abund_moorings$Mooring_presence, 
                                                  levels=unique(c("core-EGC","Shared","MIZ")))

### Plot Supplementary Figure 1a - Frequency of detection vs ASV relative abundance

Supp_Figure_1a <- ggplot(asv_max_abund_moorings) + 
  geom_point(aes(x=Num_samples_present, y=Max_rel_abund, color=Distribution_group), size=2) + 
  labs(y = "ASV maximum relative abundance (%)", x = "Number of samples ASV is present in", 
       color = "Distribution group") + 
  geom_smooth(aes(x=Num_samples_present, y=Max_rel_abund), method=lm) + 
  scale_color_manual(values=c("Intermittent" = "#C5C6D0",
                              "Transient" = "#232023", 
                              "Resident" = "#787276")) + 
  scale_y_log10() + 
  facet_wrap(Mooring_presence~., scales="free_x") + 
  guides(color=guide_legend(override.aes = list(size=4), ncol=1)) + 
  theme_bw() + 
  theme(legend.position="right", 
        legend.title = element_text(size = 16, colour="black"), 
        legend.text = element_text(size = 14, colour="black"),
        axis.text.y = element_text(size=14, colour="Black"),
        axis.title.y = element_text(size = 16, colour = "Black"),
        axis.text.x = element_text(size=14, colour="Black"),
        axis.title.x = element_text(size = 16, colour = "Black"),
        strip.background = element_rect(fill = "white"), 
        strip.text.x = element_text(size = 16, colour = "black"))
Supp_Figure_1a

### Supplementary Figure 1b
# Import data about network connections
egc_asv_connections_data=read.table("FRAM_RAS_EGC_ASV_avg_connections_network_group.txt", header = TRUE, 
                                    sep = "\t",as.is=TRUE)

## Reformat data
colnames(egc_asv_connections_data)<- gsub("\\.","-",colnames(egc_asv_connections_data))

egc_asv_connections_data_long <- melt(data=egc_asv_connections_data, id.vars=c("Group"), variable.name="Network", value.name="Connections")

egc_asv_connections_data_long$Network <- factor(egc_asv_connections_data_long$Network, 
                                                levels=unique(c("core-EGC","MIZ")))
egc_asv_connections_data_long$Group <- factor(egc_asv_connections_data_long$Group, 
                                              levels=unique(c("Intermittent", "Transient", "Resident")))

## Plot Supplementary Figure 1b
Supp_Figure_1b <- ggplot(egc_asv_connections_data_long) + 
  geom_bar(aes(x=Group, y=Connections, fill=Group), stat="identity", position="dodge", size=2) + 
  labs(y = "Avg. number of significant connections", 
       fill = "Distribution group") + 
  scale_fill_manual(values=c("Intermittent" = "#C5C6D0",
                             "Transient" = "#232023", 
                             "Resident" = "#787276")) + 
  facet_wrap(Network~., scales="free_x") + 
  guides(fill=guide_legend(override.aes = list(size=4), ncol=1)) + 
  theme_bw() + 
  theme(legend.position="right", 
        legend.title = element_text(size = 16, colour="black"), 
        legend.text = element_text(size = 14, colour="black"),
        axis.text.y = element_text(size=14, colour="Black"),
        axis.title.y = element_text(size = 16, colour = "Black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "white"), 
        strip.text.x = element_text(size = 16, colour = "black"))
Supp_Figure_1b

### Combine the above three plots into single Figure and export
pdf(file="Supplementary_Figure_1.pdf", height=8, width=12)
(Supp_Figure_1a+Supp_Figure_1b+plot_layout(widths=c(2,1), guides="collect"))
dev.off()

