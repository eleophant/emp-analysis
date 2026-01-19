########################################
#
# EMP analysis figures
#
# Eleonore Lebeuf-Taylor
#
########################################

#### packages ####
library(tidyverse)

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
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  xlab("Host sociality") +
  ylab("Observed OTUs") +
  scale_fill_viridis_d()

p_adiv_chao = df_metadata_sub %>%
  ggplot(aes(x = basic_sociality, 
             y = adiv_chao1, 
             fill = basic_sociality)) +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  xlab("Host sociality") +
  ylab("Chao1") +
  scale_fill_viridis_d()

p_adiv_shannon = df_metadata_sub %>%
  ggplot(aes(x = basic_sociality, 
             y = adiv_shannon, 
             fill = basic_sociality)) +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
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
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  xlab("Host sociality") +
  ylab("Faith PD") +
  scale_fill_viridis_d()

# combine plots
p_combined_adiv = (p_adiv_obs | p_adiv_chao) / (p_adiv_shannon | p_adiv_faith)
p_combined_adiv

theme_classic(base_size = 18)

p_combined_adiv <- p_combined_adiv +
  theme_classic(base_size = 20) +
  theme(
    strip.text = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18)
  )

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
