########################################
#
# EMP analysis figures
#
# Eleonore Lebeuf-Taylor
#
########################################

################ PACKAGES ################

library(tidyverse)
library(conflicted)
library(viridis)
library(phyloseq)
library(fantaxtic)
library(rotl)
library(phytools)
library(ape)
library(phyr)
library(ggpubr)
library(patchwork)
library(ggordiplots)
library(vegan)
library(lme4)
library(ALDEx2)
library(ggsignif)

################ DATA INPUT ################

load(file = "emp_analysis_data_2025_02_03.RData")

# run emp_analysis_revision SUBSET section to generate data for figures

################ FIGURES ################

#### Fig 1  ####  
# relative abundances
# all hosts
# normalize number of reads using median sequencing depth
total = median(sample_sums(ps_sub))
standf = function(x, t = total) round(t * (x / sum(x)))
ps_norm = transform_sample_counts(ps_sub, standf)

# mammals
total_mammals = median(sample_sums(ps_mammals_sub))
standf = function(x, t = total_mammals) round(t * (x / sum(x)))
ps_norm_mammals = transform_sample_counts(ps_mammals_sub, standf)

# relative abundances plot - mammals
# NB: to run on server only
load("plot_relab_mammals.RData")

relab_mammals_df = relab_mammals_df %>% 
  rename(Phylum = phylum)

plot_relab_mammals = ggplot(relab_mammals_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Sample") +
  ylab("Relative abundance") +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 14)) +
  scale_fill_discrete(labels = c("Actinobacteria", "Bacteroidetes", "Cyanobacteria", "Euryarchaeota", "Firmicutes", "Lentisphaerae", "Proteobacteria", "Spirochaetae",  "Tenericutes", "Verrucomicrobia", "Other")) +
  ggtitle("A (mammals)")

plot_relab_mammals

# birds
total_birds = median(sample_sums(ps_birds_sub))
standf = function(x, t = total_birds) round(t * (x / sum(x)))
ps_norm_birds = transform_sample_counts(ps_birds_sub, standf)

# convert ps to a tidy data frame
relab_birds_df = ps_norm_birds %>%
  psmelt() %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance))

# get top 10 phyla
top_phyla_birds = relab_birds_df %>%
  group_by(phylum) %>%
  summarize(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 10) %>%
  pull(phylum)

# reclassify other phyla as "Other"
relab_birds_df = relab_birds_df %>%
  mutate(phylum = ifelse(phylum %in% top_phyla_birds, phylum, "Other"))

# recalculate relative abundances with new categories
relab_birds_df = relab_birds_df %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(RelativeAbundance = Abundance / sum(Abundance)) %>% 
  rename(Phylum = phylum)

# relative abundances plot - birds
# can run on local machine
plot_relab_birds = ggplot(relab_birds_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("Sample") +
  ylab("Relative abundance") +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 14)) +
  scale_fill_discrete(labels = c("Actinobacteria", "Bacteroidetes", "Cyanobacteria", "Firmicutes", "Lentisphaerae", "Planctomycetes", "Proteobacteria", "Spirochaetae", "Tenericutes", "Verrucomicrobia", "Other")) +
  ggtitle("B (birds)")

plot_relab_birds

# Fig 1 combined
plot_relab_mammals / plot_relab_birds

#### Fig 2 ####
sociality_order = c("solitary", "intermediate", "social")
df_metadata_sub$basic_sociality = factor(df_metadata_sub$basic_sociality, levels = sociality_order)

p_adiv_obs = df_metadata_sub %>%
  ggplot(aes(x = basic_sociality, 
             y = adiv_observed_otus, 
             fill = basic_sociality)) +
  geom_boxplot() +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  xlab("Host sociality") +
  ylab("Observed OTUs") +
  scale_fill_viridis_d()

p_adiv_chao = df_metadata_sub %>%
  ggplot(aes(x = basic_sociality, 
             y = adiv_chao1, 
             fill = basic_sociality)) +
  geom_boxplot() +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  xlab("Host sociality") +
  ylab("Chao1") +
  scale_fill_viridis_d()

p_adiv_shannon = df_metadata_sub %>%
  ggplot(aes(x = basic_sociality, 
             y = adiv_shannon, 
             fill = basic_sociality)) +
  geom_boxplot() +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20)) +
  xlab("Host sociality") +
  ylab("Shannon") +
  scale_fill_viridis_d() +
  geom_signif(comparisons = list(c("social", "solitary")), 
              map_signif_level=TRUE,
              annotations="*")

p_adiv_faith = df_metadata_sub %>%
  ggplot(aes(x = basic_sociality, 
             y = adiv_faith_pd, 
             fill = basic_sociality)) +
  geom_boxplot() +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 20)) +
  xlab("Host sociality") +
  ylab("Faith PD") +
  scale_fill_viridis_d()

# combine plots
p_combined_adiv = (p_adiv_obs | p_adiv_chao) / (p_adiv_shannon | p_adiv_faith)
p_combined_adiv

ggsave(
  filename = "Figure2_300dpi.tiff",
  plot = p_combined_adiv,
  width = 330,
  height = 315,
  units = "mm",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

#### Fig 3 ####

# run emp_host_phylogeny.R

# calculate distances
dist_mat <- vegdist(df_otus_sub, method = "robust.aitchison")

# generate dispersion per species
disp_test_species <- betadisper(dist_mat, df_metadata_sub$host_species)
anova(disp_test_species) # p < 0.001

boxplot(disp_test_species, main = "Distance to centroid by social behaviour")
# create df with betadisp metrics for each species + sociality level
disp_data <- data.frame(
  distance = disp_test_species$distances,
  host_scientific_name = df_metadata_sub$host_scientific_name,
  host_species = df_metadata_sub$host_species,
  sociality = df_metadata_sub$basic_sociality
)

# filter to species with at least 5 samples
disp_data_n = disp_data %>% 
  group_by(host_scientific_name) %>% 
  summarise(n = n())

disp_data = disp_data %>% 
  left_join(disp_data_n, by = "host_scientific_name") %>% 
  filter(n > 4)

disp_data %>% 
  distinct(host_species) # 15 species, 308 observations

# add PhyPC1:5 from metadata with n > 4
df_metadata_phy = df_metadata_sub %>% 
  select(host_species, host_scientific_name, host_genus, starts_with("Phy")) %>% 
  distinct(host_scientific_name, .keep_all = TRUE) # keep all columns

disp_data = disp_data %>%
  left_join(df_metadata_phy, by = "host_scientific_name")

# set reference level to 'social'
disp_data$sociality <- factor(disp_data$sociality, levels = c("social", "intermediate", "solitary"))

# test mean dispersion per species, weighted by sample size bc more samples = more reliable estimate
species_disp <- disp_data %>%
  group_by(host_scientific_name, sociality) %>%
  summarise(mean_distance = mean(distance), n = n())

lm_weighted <- lm(mean_distance ~ sociality, data = species_disp, weights = n)
summary(lm_weighted)

# ordering
disp_data = disp_data %>% 
  mutate(sociality = fct_relevel(sociality, c("solitary", "intermediate", "social"))) %>%
  arrange(sociality, n) %>%  # reorder by sociality then sample size
  mutate(host_scientific_name = factor(host_scientific_name, levels = unique(host_scientific_name)))

# plot dispersion distances
fig3 <- disp_data %>%
  ggplot(aes(x = host_scientific_name, y = distance, fill = sociality)) +
  geom_boxplot(width = 0.6) +
  scale_fill_viridis_d(name = "Sociality", direction = -1) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 140, by = 20)) +
  geom_text(data = disp_data, aes(host_scientific_name, Inf, label = n), hjust = "inward") +
  labs(x = "Host species", y = "Distance from centroid") +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip()

ggsave(
  filename = "Figure3_300dpi.tiff",
  plot = fig3,
  width = 330,
  height = 315,
  units = "mm",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)

#### Fig 4 ####
# run BETA DIV code snippet

plot_rda_all = plot_rda_all + theme(legend.position = "none")
plot_disp_sociality = plot_disp_sociality + theme(legend.position = "none")
plot_disp_soc_boxplot + theme(legend.position = "none")
fig4AB = plot_rda_all + plot_disp_sociality # A + B
fig4C = plot_disp_soc_boxplot # C
fig4_combined = fig4AB / fig4C

library(gridExtra)

# convert plots to grobs
g1 <- ggplotGrob(plot_rda_all)
g2 <- ggplotGrob(plot_disp_sociality)
g3 <- ggplotGrob(plot_disp_soc_boxplot)

# combine A + B horizontally
top_row <- arrangeGrob(g1, g2, ncol = 2)
# combine top row and C vertically
fig4_combined <- arrangeGrob(top_row, g3, ncol = 1)

# save
ggsave(
  filename = "Figure4_300dpi.tiff",
  plot = fig4_combined,
  width = 300, 
  height = 300, 
  units = "mm",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)


#### Fig 5 ####

ggsave(
  filename = "Figure5_300dpi.tiff",
  plot = plot_prosocials,
  width = 300,
  height = 200,
  units = "mm",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)
