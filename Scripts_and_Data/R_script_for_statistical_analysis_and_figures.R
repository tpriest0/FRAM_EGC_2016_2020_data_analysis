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

setwd("D:/ownCloud/PhD/Arctic Bacterioplankton/Fram Strait/RAS_metagenomes/EGC/ASV/FRAM_RAS_EGC_data_for_processing")


############################################################################################
### LOADING PACKAGES AND DATA ###
############################################################################################

BiocManager::install("ALDEx2", lib="D:/R")
options(repos = c(
  fawda123 = 'https://fawda123.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

install.packages('ggord')

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
library(ggpmisc)

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
  geom_point(aes(x=Mooring_name,y=value), colour="black", size = 2) + 
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

# Import network connection information
egc_asv_network_connections=read.table("FRAM_RAS_EGC_ASV_avg_connections_network_group.txt", header = TRUE, 
                        sep = "\t",check.names=F,as.is=TRUE)

###
# Plot Figure 3a (asv dynamics)

egc_asv_dynamics_df$Mooring_presence <- factor(egc_asv_dynamics_df$Mooring_presence, 
                                                 levels=unique(c("core-EGC", "Shared","MIZ")))

Figure_3a <- ggplot(egc_asv_dynamics_df, aes(x=Num_samples_present, y=Max_rel_abund)) + 
  geom_point(aes(colour=Distribution_group), stat="identity", show.legend=FALSE) + 
  stat_poly_eq(method="lm", use_label(c("R2", "P")),
               label.y=0.95) + 
  stat_poly_line() + 
  labs(y = "ASV maximum relative abundance (%)", x = "Number of samples ASV is present in",
       colour = "Distribution group") + 
  scale_y_log10() + 
  facet_grid(.~Mooring_presence, scales="free_x") + 
  scale_colour_manual(values=c("Intermittent" = "#C5C6D0",
                             "Transient" = "#232023", 
                             "Resident" = "#787276")) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12, colour="black"),
        axis.title.y = element_text(size=14, colour="black"),
        axis.text.x = element_text(size=12, colour="black"),
        axis.title.x = element_text(size=14, colour="black"),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 14, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        strip.background = element_rect(fill="white"), 
        strip.text.x = element_text(size = 14, colour = "black"))
Figure_3a

### 
# Reformat table on network connections
egc_asv_network_connections_long = egc_asv_network_connections %>%
  reshape2::melt(., id.vars="Group", variable.name="Mooring", 
                 value.name="num_connections")
egc_asv_network_connections_long$Mooring <- factor(egc_asv_network_connections_long$Mooring, 
                                               levels=unique(c("core-EGC","MIZ")))
egc_asv_network_connections_long$Group <- factor(egc_asv_network_connections_long$Group, 
                                                   levels=unique(c("Intermittent","Transient",
                                                                   "Resident")))

# Plot Figure 3b (network connections per group)

Figure_3b <- ggplot(egc_asv_network_connections_long, aes(x = Group, y = num_connections)) + 
  geom_bar(aes(fill=Group), stat="identity", position="identity") + 
  labs(y = "Avg. number of significant connections", fill = "Distribution group") +
  facet_grid(.~Mooring, scales="free") + 
  scale_fill_manual(values=c("Intermittent" = "#C5C6D0",
                               "Transient" = "#232023", 
                               "Resident" = "#787276")) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12, colour="black"),
        axis.title.y = element_text(size=14, colour="black"),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size=14, colour="black"),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 14, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        strip.background = element_rect(fill="white"), 
        strip.text.x = element_text(size = 14, colour = "black"))

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

### 
# Plot Figure 3c - Dynamics of community fractions across samples

Figure_3c <- ggplot(data=distr_group_rel_abund, 
                                          aes(x=as.Date(Date), y=Rel_abund)) + 
  scale_x_date(breaks = "4 months") + 
  labs(x = "Date", y = "Relative abundance (%)", fill = "Distribution group") + 
  geom_bar(aes(fill=Distribution_group), stat="identity", position="stack", width = 8) + 
  facet_grid(.~Mooring_position, scales="free_x") + 
  scale_fill_manual(values=c("Intermittent" = "#C5C6D0",
                             "Transient" = "#232023", 
                             "Resident" = "#787276")) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size=12, colour="Black"),
        axis.title.y = element_text(size=14, colour="Black"),
        axis.text.x = element_text(size=12, colour="Black", angle=45, hjust=1),
        axis.title.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_text(size = 14, colour = "black"), 
        legend.text = element_text(size = 12, colour = "black"), 
        strip.background = element_rect(fill="white"), 
        strip.text.x = element_text(size = 14, colour = "black"))
Figure_3c

pdf(file = "Figure_3.pdf", height=8, width=12)
((Figure_3a+Figure_3b)+plot_layout(widths=c(2,1)))/
  Figure_3c
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

### Figure 5. Comparison and assessment of distribution of MAGs across Arctic Ocean
### and Fram Strait metagenomes 

##########

# Load libraries
library(dplyr)
library(tidyverse)
library(reshape2)
library(vegan)
library(UpSetR)

### Part A - Species-level MAG comparison

# Set working directory
setwd('D:/ownCloud/PhD/Arctic Bacterioplankton/Arctic_comparison/MAGs/Arctic_Fram_EGC/')

# Load MAG species clusters summary info
comparative_mag_table <- read.table(
  "FRAM_RAS_EGC_vs_Arctic_MAGs_species_cluster_summary.txt",
  h = T, sep = "\t",
  check.names=F)


# Reformat dataframe ready for upset plot
MAG_species_cluster_info_for_upset = comparative_mag_table %>%
  subset(., select=c(Data_source,Species)) %>%
  mutate(count=1) %>%
  reshape2::dcast(., Species~Data_source, value.var="count")%>%
  mutate_if(is.numeric, ~1 * (. >0)) %>%
  subset(., select=-c(Species))

datasets_list = colnames(MAG_species_cluster_info_for_upset)

temp1 <- as.data.frame(datasets_list) %>%
  rename(sets = datasets_list)
temp2 <- as.data.frame(datasets_list) %>%
  rename(dataset_source = datasets_list)

upset_metadata <- cbind(temp1,temp2)

# Plot Figure 5a
Figure_5a <- 
  UpSetR::upset(MAG_species_cluster_info_for_upset, order.by = "freq", decreasing = T,
                sets.x.label = "Number of species", 
                mainbar.y.label = "Number of species",
                text.scale = c(1.5,1,1.5,1,1.5,1),
                set.metadata = list(data=upset_metadata, 
                                    plots = list(list(type = "matrix_rows",
                                                      column = "dataset_source",
                                                      colors=c(FRAM_EGC="#4D4D5C",
                                                               FRAM18="#C5C6D0",
                                                               MOSAiC="#4D4D5C",
                                                               TARA="#C5C6D0",
                                                               alpha=0.5)))))


#####

### Part B - Frequency of detection of FRAM_EGC MAGs across metagenomes 

# Set working directory
setwd('D:/ownCloud/PhD/Arctic Bacterioplankton/Arctic_comparison/MAGs/Arctic_Fram_EGC/')

# Load MAG species relative abundance file
mag_species_abund_arctic_fram_long <- read.table(
  "Arctic_Fram_MAGs_species_reps_rel_abund_long.txt",
  h = T, sep = "\t",
  check.names=F)


# Reformat
mag_species_abund_arctic_fram_long %>%
  subset(., select="Sample") %>%
  count(Sample)

# Load MAG species clusters summary info
egc_mag_species_info <- read.table(
  "Arctic_Fram_MAGs_EGC_reps_cluster_info.txt",
  h = T, sep = "\t",
  check.names=F)

# Extract only species for which an FRAM_EGC MAG was recovered and reformat
# dataframe to wide
mag_species_abund_arctic_fram_wide = mag_species_abund_arctic_fram_long %>%
  filter(Species %in% egc_mag_species_info$Species) %>%
  subset(., select=c(Sample,Species,Rel_abund)) %>%
  aggregate(Rel_abund~Species+Sample, data=., FUN=sum) %>%
  reshape2::dcast(., Species~Sample, value.var = "Rel_abund") %>%
  tibble::column_to_rownames(., var="Species")

# Extract only species for which an FRAM_EGC MAG was recovered and reformat
# dataframe to wide and calculate frequency of detection across all samples
species_freq_detection = mag_species_abund_arctic_fram_long %>%
  filter(Species %in% egc_mag_species_info$Species) %>%
  subset(., select=c(Sample,Species,TAD)) %>% 
  aggregate(TAD~Species+Sample, data=., FUN=sum) %>%
  reshape2::dcast(., Species~Sample, value.var = "TAD") %>%
  tibble::column_to_rownames(., var="Species")%>%
  replace(., . > 1, 1) %>% 
  replace(., . < 1, 0) %>%
  mutate(Num_samples_present = rowSums(.[1:67])) %>% 
  as.data.frame(.) %>% 
  tibble::rownames_to_column(., var="Species") %>%
  subset(., select=c(Species,Num_samples_present))

# Determine the min, mean and max relative abundance value for each MAG
# MAX
MAG_max_rel_abund = mag_species_abund_arctic_fram_wide %>% 
  tibble::rownames_to_column(., var="Species") %>%
  melt(., id.vars="Species", variable.name="Sample", value.name="Rel_abund") %>% 
  aggregate(Rel_abund~Species, data=., max) %>% 
  rename(., Max_rel_abund = Rel_abund)
# MIN
MAG_min_rel_abund = mag_species_abund_arctic_fram_wide %>% 
  tibble::rownames_to_column(., var="Species") %>%
  melt(., id.vars="Species", variable.name="Sample", value.name="Rel_abund") %>% 
  aggregate(Rel_abund~Species, data=., min) %>% 
  rename(., Min_rel_abund = Rel_abund)
# MEAN
MAG_avg_rel_abund = mag_species_abund_arctic_fram_wide %>% 
  tibble::rownames_to_column(., var="Species") %>%
  melt(., id.vars="Species", variable.name="Sample", value.name="Rel_abund") %>% 
  aggregate(Rel_abund~Species, data=., mean) %>% 
  rename(., Avg_rel_abund = Rel_abund)

# Combine different statistics generated above into a single summary tabl
egc_mag_arctic_summary_info = egc_mag_species_info %>%
  left_join(species_freq_detection, by="Species") %>%
  left_join(MAG_max_rel_abund, by="Species") %>%
  left_join(MAG_avg_rel_abund, by="Species") %>%
  left_join(MAG_min_rel_abund, by="Species") %>%
  arrange(desc(Num_samples_present))

# Combine with summary information on MAG and its connecting ASV generated earlier 
# in this script
setwd('D:/ownCloud/PhD/Arctic Bacterioplankton/Fram Strait/RAS_metagenomes/EGC/ASV/FRAM_RAS_EGC_data_for_processing')

# Import files if not running script consecutively
mag_to_asv_mapping=read.table("FRAM_RAS_EGC_ASV_to_MAG.txt", header = TRUE, 
                              sep = "\t",as.is=TRUE,check.names=F)
asv_summary_complete=read.table("FRAM_RAS_EGC_ASV_dynamics_summary_and_taxa_complete_with_spls_cluster.txt", header = TRUE, 
                                sep = "\t",as.is=TRUE,check.names=F)

setwd('D:/ownCloud/PhD/Arctic Bacterioplankton/Fram Strait/RAS_metagenomes/EGC/MAGs/')

mag_reps_summary=read.table("FRAM_RAS_EGC1617_MAG_reps_summary.txt", header = TRUE, 
                            sep = "\t",as.is=TRUE,check.names=F)

# Setwd ready for plot output
setwd('D:/ownCloud/PhD/Arctic Bacterioplankton/manuscripts/FRAM_RAS_EGC1620_manuscript/ISME_submission/Resubmission_after_review/working_figures')

# Reformat naming of samples to match across dataframes
mag_reps_summary_mod = mag_reps_summary %>%
  mutate(across("MAG", str_replace, "_RAS_EGC1617", "_EGC1617"))

# Extract the information of interest from asv summary dataframe
asv_distribution_group_and_cluster = asv_summary_complete %>%
  subset(., select=c(ASV_name,Distribution_group,cluster))

# Combine MAG distribution statistics with taxonomic information and
# linked ASV 
egc_mag_summary_across_arctic_and_fram = egc_mag_arctic_summary_info %>%
  left_join(mag_reps_summary_mod, by="MAG") %>%
  left_join(mag_to_asv_mapping, by=c("MAG" = "MAG_name")) %>% 
  left_join(asv_distribution_group_and_cluster, by="ASV_name") %>%
  subset(., select=c(MAG,Species,Rank,Num_samples_present,Completeness,
                     Contamination,Strain_heterogeneity,Genome_size,
                     N50,GC_content,GTDB_taxonomy,MAG_sample_source,ASV_name,
                     Distribution_group,Max_rel_abund,Avg_rel_abund,
                     Min_rel_abund,cluster)) %>%
  mutate(Distribution_group = Distribution_group %>%
           is.na %>% 
           ifelse("Unassigned",Distribution_group)) %>% 
  mutate(Distribution_group = fct_relevel(Distribution_group, c("Resident",
                                                                "Intermittent",
                                                                "Transient",
                                                                "Unassigned"))) %>%
  arrange(Distribution_group,desc(Num_samples_present))


# Export MAG dynamics summary information file
write.table(egc_mag_summary_across_arctic_and_fram,
            file="FRAM_RAS_EGC_MAGs_dynamics_summary_information.txt", sep="\t")

# Check for statistical significance in the frequency of detection between
# different distribution groups.
# Need to remove Transient however, due to it only having one MAG rep
distr_group_num_samples = egc_mag_summary_across_arctic_and_fram %>%
  subset(., select=c(Num_samples_present,Distribution_group)) %>%
  filter(Distribution_group != "Transient")
head(distr_group_num_samples)

kruskal.test(Num_samples_present~Distribution_group, data=distr_group_num_samples) 
pairwise.wilcox.test(distr_group_num_samples$Num_samples_present, 
                     distr_group_num_samples$Distribution_group, p.adjust.method = "BH")

# Fix MAG order
egc_mag_summary_across_arctic_and_fram$MAG <- 
  factor(egc_mag_summary_across_arctic_and_fram$MAG,
         levels=unique(egc_mag_summary_across_arctic_and_fram$MAG))

# Plot Figure 5b
Figure_5b = 
  ggplot(data=egc_mag_summary_across_arctic_and_fram,
         aes(x = Distribution_group, y = Num_samples_present)) + 
  geom_boxplot(aes(fill = Distribution_group)) + 
  geom_jitter(alpha = 0.8, width = 0.2) +
  geom_hline(aes(yintercept=67), colour="purple", linewidth=1.5) + 
  labs(y = "Number of samples MAG has >1X coverage") + 
  scale_fill_manual(values=c("Resident" = "#787276",
                             "Intermittent" = "#C5C6D0",
                             "Transient" = "#232023",
                             "Unassigned" = "#747D9C")) +
  theme_bw() + 
  theme(legend.position="none", 
        axis.title.y = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(colour = "black", size = 12))
Figure_5b

#####

### Part C - Dynamics of sPLS clusters across Fram Strait and Arctic Ocean 
### metagenomes

# Set working directory
setwd('D:/ownCloud/PhD/Arctic Bacterioplankton/Arctic_comparison/')

fram_arctic_mg_metadata <- read.table(
  "FRAM_Arctic_MG_all_metadata.txt",
  h = T, sep = "\t",
  check.names=F)

# Extract necessary information from mag summary file generated above
egc_mag_name_to_species = egc_mag_summary_across_arctic_and_fram %>%
  subset(., select=c(MAG,Species,Distribution_group,Num_samples_present))
View(egc_mag_name_to_species)

egc_mag_to_cluster = subset(egc_mag_summary_across_arctic_and_fram, select=c(MAG,cluster))

# Create average relative abundance dataframe for spls clusters across samples and
# combine with metadata
spls_clusters_abund_across_arctic_long_avg = mag_species_abund_arctic_fram_long %>% 
  #filter(TAD > 1) %>% 
  left_join(egc_mag_name_to_species, by="Species") %>%
  left_join(egc_mag_to_cluster, by="MAG") %>%
  aggregate(Rel_abund~cluster+Sample, data=., FUN=mean) %>%
  left_join(fram_arctic_mg_metadata, by = "Sample") %>%
  arrange(desc(Depth)) %>%
  mutate(Dataset = case_when(
    grepl("TARA",Sample)~"TARA",
    grepl("MOSAIC", Sample)~"MOSAIC",
    grepl("EGC",Sample)~"FRAM_EGC",
    grepl("FRAM18",Sample)~"FRAM18"
  ))

# Define colour scheme for clusters
pls_cluster_colours <- c(
  "C1" = "#db6d00",
  "C2" = "#920000",
  "C3" = "#b6dbff",
  "C4" = "#006ddb",
  "C5" = "#009292",
  "C6" = "#ffb6db",
  "C7" = "#ffff6d",
  "C8" = "#000000")

### Plot Figure 5c
Figure_5c_ice <- ggplot(spls_clusters_abund_across_arctic_long_avg, 
                        aes(y=Rel_abund, x=Ice_cover, fill=cluster)) +
  #geom_area(aes(y=Ice_cover, x=ifelse(Rel_abund>0.05,Rel_abund,NA), fill=cluster)) + 
  scale_fill_manual(values=pls_cluster_colours) + 
  stat_smooth(se = FALSE, geom = 'area', method = 'loess', span = 1, alpha=0.5, 
              position = "stack", aes(fill = cluster)) + 
  scale_y_continuous(limits = c(0,3)) + 
  labs(y = "Average relative abundance (%)", x = "Ice cover (%)") + 
  theme_bw() + 
  theme(legend.position = "right", 
        legend.title = element_text(colour = "black", size = 14),
        legend.text = element_text(colour = "black", size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.text.x = element_text(colour = "black", size = 12))
Figure_5c_ice

# Plot spls abundance area with daylight
Figure_5c_daylight <- ggplot(data=spls_clusters_abund_across_arctic_long_avg, 
                             aes(y=ifelse(Rel_abund>0,Rel_abund,NA), x=Daylight, fill=cluster)) +
  #geom_area(aes(y=Ice_cover, x=ifelse(Rel_abund>0.05,Rel_abund,NA), fill=cluster)) + 
  scale_fill_manual(values=pls_cluster_colours) + 
  stat_smooth(se = FALSE, geom = 'area', method = 'loess', span = 1, alpha=0.5, 
              position = "stack", aes(fill = cluster)) + 
  scale_y_continuous(limits = c(0,4.5)) + 
  labs(y = "Average relative abundance (%)", x = "Daylight (h)") + 
  theme_bw() + 
  theme(legend.position = "right", 
        legend.title = element_text(colour = "black", size = 14),
        legend.text = element_text(colour = "black", size = 12), 
        axis.title.y = element_blank(),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.text.x = element_text(colour = "black", size = 12))
Figure_5c_daylight

# Plot spls abundance area with depth
Figure_5c_depth = spls_clusters_abund_across_arctic_long_avg %>%
  filter(Depth < 1000) %>%
  ggplot(data=., aes(y=Rel_abund, x=Depth, fill=cluster)) +
  stat_smooth(se = FALSE, geom = 'area', method = 'loess', span = 1, alpha=0.5, 
              position = "stack", aes(fill = cluster)) + 
  scale_fill_manual(values=pls_cluster_colours) + 
  scale_y_continuous(limits = c(0,4)) + 
  labs(y = "Average relative abundance (%)", x = "Depth (m)") +
  theme_bw() + 
  theme(legend.position = "right", 
        legend.title = element_text(colour = "black", size = 14),
        legend.text = element_text(colour = "black", size = 12), 
        axis.title.y = element_text(colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 12),
        axis.title.x = element_text(colour = "black", size = 14),
        axis.text.x = element_text(colour = "black", size = 12))
Figure_5c_depth

### Combine all Figure 5 parts into one
library(cowplot)
library(ggplotify)

fram_arctic_MAG_species_comparison_upset_v2 <- as.grob(fram_arctic_MAG_species_comparison_upset)

part_1 <- plot_grid(fram_arctic_MAG_species_comparison_upset_v2, 
                    FRAM_EGC_MAGs_distribution_group_across_arctic_boxplot, 
                    nrow = 2, align="hv", labels = c("A", "B"),
                    label_size = 14)

part_2 <- plot_grid(spls_clust_avg_abund_with_ice,
                    spls_clust_avg_abund_with_daylight,
                    spls_clust_avg_abund_with_depth,
                    nrow = 3, align = "v", labels = "C", label_size = 14)

# Setwd ready for plot output
setwd('D:/ownCloud/PhD/Arctic Bacterioplankton/manuscripts/FRAM_RAS_EGC1620_manuscript/ISME_submission/Resubmission_after_review/working_figures')

pdf(file="Figure_5.pdf", height = 10, width = 10)
plot_grid(part_1, part_2, ncol=2)
dev.off()



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
pdf(file="Figure_6.pdf", height=12, width=14)
(egc_asv_meta_summary_plot/pls_c1_temporal/pls_c2_temporal/
  pls_c3_temporal/pls_c4_temporal/pls_c6_temporal/pls_c7_temporal/pls_c8_temporal)+
  plot_layout(guides="collect")
dev.off()

### Plot was further modified in Inkscape


############################################################################################
### Figure 7
############################################################################################

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
### Gene-encoding functions were checked and rename/grouped based on function


##########

### Figure 8. Selected genes involved in the uptake and degradation 
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
    fill=cluster), outlier.shape=NA) + 
  geom_point(aes(
    x=Gene_short, y=CLR_value, colour=cluster), position=position_dodge(width=0.75)) +
  geom_hline(
    yintercept=0, size=0.5,
    linetype="dashed", 
    color="gray33") + 
  scale_fill_manual(values=c(
    "High-ice cluster"="darkcyan",
    "Low-ice cluster"="lightgoldenrod2")) +
  scale_colour_manual(values=c(
    "High-ice cluster"="black",
    "Low-ice cluster"="black")) + 
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

pdf(file="Figure_7.pdf", height = 12, width = 12)
aldex_enriched_genes_expanded_boxplot
dev.off()

############################################















