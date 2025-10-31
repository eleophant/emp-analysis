########################################
#
# EMP analysis on subsampled data
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
theme_set(theme_classic())
filter = dplyr::filter
select = dplyr::select
conflicts_prefer(dplyr::slice)

################ DATA INPUT ################

load(file = "emp_analysis_data_2025_02_03.RData") # this dataset has been modified to omit 1 reptile species, black bear, colobine primates

################ SUBSET ################

#### make dfs ####
df_obs_n = df_metadata %>% 
  group_by(host_species) %>% 
  summarise(n = n()) %>% 
  print(n = 27)

df_obs_n %>% summarise(mean(n), median(n), max(n), min(n))
# mean 41, median 3, max 902, min 1

# subsample metadata
df_metadata_sub = df_metadata %>% 
  rownames_to_column() %>% 
  mutate(sample_id = rowname) %>%  # keep a sample_id column
  group_by(host_species) %>% 
  slice_sample(n = 41)  %>%
  ungroup() %>% 
  column_to_rownames("rowname")

# subsample otu table
common_samples_sub = intersect(rownames(mx_otus), rownames(df_metadata_sub))

df_otus_sub = mx_otus %>% as.data.frame() %>% 
  rownames_to_column() %>% 
  filter(rowname %in% common_samples_sub) %>% 
  column_to_rownames("rowname") %>%
  select(where(~ sum(.) != 0)) # remove empty columns (OTUs that are 0 everywhere)

# ensure that row order matches
df_metadata_sub = df_metadata_sub %>%
  slice(match(rownames(df_otus_sub), rownames(df_metadata_sub)))
all(rownames(df_otus_sub) == rownames(df_metadata_sub)) # TRUE

#### phyloseq ####
# set up metadata & otu table
p_metadata_sub = sample_data(df_metadata_sub)
p_otu_sub = otu_table(df_otus_sub, taxa_are_rows = FALSE)

# check that sample names match
all(sample_names(p_metadata_sub) == sample_names(p_otu_sub)) # TRUE

# assemble phyloseq object
physeq_sub = phyloseq(p_otu_sub, p_metadata_sub, p_tax)
ps_sub = merge_phyloseq(physeq_sub, tree_emp)

#### subset classes ####
# mammals
df_metadata_sub_mammals <- df_metadata_sub |>
  filter(host_class == "c__Mammalia")

# keep only samples in OTU matrix
common_samples_sub_mammals <- intersect(rownames(mx_otus), rownames(df_metadata_sub_mammals))

# order-preserving approach
df_metadata_sub_mammals <- df_metadata_sub_mammals[common_samples_sub_mammals, , drop = FALSE]

df_otus_sub_mammals <- mx_otus |>
  as.data.frame() |>
  rownames_to_column(var = "rowname") |>
  filter(rowname %in% common_samples_sub_mammals) |>
  column_to_rownames("rowname") |>
  select(where(~ sum(.) != 0))

# force both to the same order
df_otus_sub_mammals <- df_otus_sub_mammals[rownames(df_metadata_sub_mammals), ]

# check alignment
all(rownames(df_otus_sub_mammals) == rownames(df_metadata_sub_mammals)) # TRUE

# birds
df_metadata_sub_birds = df_metadata_sub |>
  filter(host_class == "c__Aves")

common_samples_sub_birds = intersect(rownames(mx_otus), rownames(df_metadata_sub_birds))

df_otus_sub_birds = mx_otus |> as.data.frame() |> 
  rownames_to_column() |> 
  filter(rowname %in% common_samples_sub_birds) |> 
  column_to_rownames("rowname") |>
  select(where(~ sum(.) != 0))

all(rownames(df_otus_sub_birds) == rownames(df_metadata_sub_birds)) # TRUE

# phyloseq
ps_mammals_sub = subset_samples(ps_sub, host_class == "c__Mammalia")
ps_birds_sub = subset_samples(ps_sub, host_class == "c__Aves")

################ SEQ DEPTH ################

#### setup ####
sequencing_depth <- rowSums(df_otus_sub)

cat("=== Sequencing Depth Summary ===\n")
summary(sequencing_depth)

cat("\nStandard deviation:", sd(sequencing_depth), "\n")
cat("Coefficient of variation:", sd(sequencing_depth) / mean(sequencing_depth), "\n")

cat("\nRange:", min(sequencing_depth), "to", max(sequencing_depth), "\n")
cat("Fold difference:", max(sequencing_depth) / min(sequencing_depth), "\n")

# create dataframe with depth
depth_data = data.frame(
  sample_id = rownames(df_otus_sub),
  depth = rowSums(df_otus_sub)
) %>%
  left_join(
    df_metadata_sub %>% 
      select(sample_id, host_scientific_name, basic_sociality, study_id),
    by = "sample_id"
  )

# plot histogram
ggplot(depth_data, aes(x = depth)) +
  geom_histogram(bins = 50, fill = "skyblue2", color = "black") +
  geom_vline(xintercept = median(depth_data$depth), 
             color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Distribution of Sequencing Depth",
    subtitle = paste("Median =", round(median(depth_data$depth)), "reads"),
    x = "Sequencing depth (total reads)",
    y = "Number of samples"
  ) +
  theme_minimal()
# looks highly skewed (long tail)

# log scale version
ggplot(depth_data, aes(x = depth)) +
  geom_histogram(bins = 50, fill = "skyblue2", color = "black") +
  scale_x_log10() +
  labs(
    title = "Distribution of Sequencing Depth (Log Scale)",
    x = "Sequencing depth (log10)",
    y = "Number of samples"
  ) +
  theme_minimal()

### alpha div ####

# check if depth correlates with alpha diversity
depth_alpha <- data.frame(
  depth = sequencing_depth,
  adiv_observed_otus = df_metadata_sub$adiv_observed_otus,
  adiv_chao1 = df_metadata_sub$adiv_chao1,
  adiv_shannon = df_metadata_sub$adiv_shannon,
  adiv_faith_pd = df_metadata_sub$adiv_faith_pd,
  sociality = df_metadata_sub$basic_sociality
)

# correlation test
cor_test_obs_otus <- cor.test(depth_alpha$depth, depth_alpha$adiv_observed_otus, method = "spearman")

cor_test_chao1 <- cor.test(depth_alpha$depth, depth_alpha$adiv_chao1, method = "spearman")

cor_test_shannon <- cor.test(depth_alpha$depth, depth_alpha$adiv_shannon, method = "spearman")

cor_test_faith_pd <- cor.test(depth_alpha$depth, depth_alpha$adiv_faith_pd, method = "spearman")

cat("Spearman correlations between depth and alpha diversity:\n")

print(cor_test_obs_otus) # rho = 0.02309349, p = 0.6558
print(cor_test_chao1) # rho = 0.01962841, p = 0.7048
print(cor_test_shannon) # rho = 0.005850495, p = 0.9101
print(cor_test_obs_otus) # rho = 0.02309349, p = 0.6558


# plot
ggplot(depth_alpha, aes(x = depth, y = alpha_div, color = sociality)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_x_log10() +
  labs(
    title = "Alpha diversity vs Sequencing Depth",
    subtitle = paste("Spearman ρ =", round(cor_test$estimate, 3), 
                     ", p =", format.pval(cor_test$p.value, digits = 3)),
    x = "Sequencing depth (log10)",
    y = "Observed OTUs",
    color = "Sociality"
  ) +
  theme_minimal()


#### sociality ####
ggplot(depth_data, aes(x = basic_sociality, y = depth, fill = basic_sociality)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  scale_y_log10() +
  labs(
    title = "Sequencing depth by sociality level",
    x = "Sociality",
    y = "Sequencing depth (log10)",
    fill = "Sociality"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
# looks like solitary has lower sequencing depth

kruskal.test(depth ~ basic_sociality, data = depth_data)
# chi-squared = 53.15, df = 2, p = 2.875e-12

# pairwise comparisons
pairwise.wilcox.test(depth_data$depth, depth_data$basic_sociality, p.adjust.method = "BH")
# sig: solitary vs social, solitary vs intermediate


################ DATA EXPLORATION ################

### collinearity ####

table_sd <- table(df_metadata_sub$basic_sociality, df_metadata_sub$basic_diet)

chisq.test(table_sd)
# X-squared = 101.75, df = 4, p-value < 2.2e-16
# soiality and diet are highly collinear

#### all hosts ####
df_metadata_sub %>% nrow() # 375 samples
df_metadata_sub %>% group_by(host_species) %>% 
  summarise(n = n()) %>% print(n = 100) # 45 host species
df_otus_sub %>% as.data.frame() %>% ncol() # number of OTUs = 34,482
df_otus_sub %>% sum() # total reads = 13,085,267

# per sample
mean(rowSums(df_otus_sub)) # mean number of reads = 34,894.05
sd(rowSums(df_otus_sub)) # 24,286.86

otus_per_sample = df_otus_sub %>% mutate(total_otus = rowSums(df_otus_sub != 0)) %>% rownames_to_column() %>% 
  select(rowname, total_otus)

otus_per_sample = df_metadata_sub %>% 
  rownames_to_column() %>% 
  # rename(rowname = sample_id) %>% 
  full_join(otus_per_sample) %>% 
  select(rowname, host_class, host_scientific_name, total_otus)

otus_per_sample %>% summarise(mean(total_otus), median(total_otus), max(total_otus), min(total_otus))
# mean 1019, median 961, max 2817, min 114

#### mammals ####
df_metadata_sub_mammals %>% nrow() # 312 samples
df_metadata_sub_mammals %>% group_by(host_species) %>% 
  summarise(n = n()) %>% print(n = 35) # 35 host species
df_otus_sub_mammals %>% ncol() # 31,451 OTUs
df_otus_sub_mammals %>% sum() # total reads = 11,086,194

# per sample
mean(rowSums(df_otus_sub_mammals)) # mean number of reads = 35,532.67
sd(rowSums(df_otus_sub_mammals)) # 23,866.22

otus_per_sample_mammals = df_otus_sub_mammals %>% mutate(total_otus = rowSums(df_otus_sub_mammals != 0)) %>% rownames_to_column() %>% 
  select(rowname, total_otus)

otus_per_sample_mammals = df_metadata_sub_mammals %>% 
  # rename(rowname = sample_id) %>% 
  rownames_to_column() %>% 
  full_join(otus_per_sample_mammals) %>% 
  select(rowname, host_class, host_scientific_name, total_otus)

otus_per_sample_mammals %>% summarise(mean(total_otus), median(total_otus), max(total_otus), min(total_otus))
# mean 1062, median 1004, max 2473, min 209

#### birds ####
df_metadata_sub_birds %>% nrow() # 63 samples
df_metadata_sub_birds %>% group_by(host_species) %>% 
  summarise(n = n()) %>% print(n = 100) # 10 host species
df_otus_sub_birds %>% ncol() # 13,915
df_otus_sub_birds %>% sum() # total reads = 1,999,073

# per sample
mean(rowSums(df_otus_sub_birds)) # mean number of reads = 31,731.32
sd(rowSums(df_otus_sub_birds)) # 26,244.05

# per sample
otus_per_sample_birds = df_otus_sub_birds %>% mutate(total_otus = rowSums(df_otus_sub_birds != 0)) %>% rownames_to_column() %>% 
  select(rowname, total_otus)

otus_per_sample_birds = df_metadata_sub_birds %>% 
  # rename(rowname = sample_id) %>% 
  rownames_to_column() %>% 
  full_join(otus_per_sample_birds) %>% 
  select(rowname, host_class, host_scientific_name, total_otus)

otus_per_sample_birds %>% summarise(mean(total_otus), median(total_otus), max(total_otus), min(total_otus))
# mean 807, median 695, max 2817, min 114

#### compare classes ####
wilcox.test(total_otus ~ host_class, data = otus_per_sample) # W = 6841, p-value = 0.0001415

#### social OTUs ####
# question: are (Limosi)lactobacillus and Bifidobacterium present in df_otus?
# extract otu tax table from phyloseq object
otu_tax = ps_sub %>% tax_table()
df_tax = otu_tax %>%
  as(., "matrix") %>%
  as.data.frame() %>%
  rownames_to_column()

# Bifidobacterium
df_tax %>% filter(genus == "D_5__Bifidobacterium") %>% 
  summarise(n = n()) # 14 Bifido species

# make list of Bifidobacterium OTUs
ls_bifido = df_tax %>%
  filter(genus == "D_5__Bifidobacterium") %>% 
  pull(rowname)

#Limosilactibacillus
df_tax %>% filter(genus == "D_5__Lactobacillus") %>% 
  summarise(n = n()) # 226 Limosilactobacillus species

# make list of Limosilactobacillus OTUs
ls_limosi = df_tax %>%
  filter(genus == "D_5__Lactobacillus") %>% 
  pull(rowname)

# total OTUs
df_otus_sub %>% 
  sum() %>% # total = 13,085,267

# Bifidobact OTUs
df_otus_sub %>% 
  select(all_of(ls_bifido)) %>%
  summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>% 
  sum() # 2685
# Bifido: 2685/13085267*100 = 0.02051926 % of sequences

# Limosilact OTUs
df_otus_sub %>% 
  select(all_of(ls_limosi)) %>% # all Limosilact OTUs
  summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>% 
  sum() # 59335
# Limosilact: 59335/13085267*100 = 0.4534489 % of sequences

#### relative abundances ####  
####_ all hosts #### 
# normalize number of reads using median sequencing depth
total = median(sample_sums(ps_sub))
standf = function(x, t = total) round(t * (x / sum(x)))
ps_norm = transform_sample_counts(ps_sub, standf)

####  _ mammals #### 
total_mammals = median(sample_sums(ps_mammals_sub))
standf = function(x, t = total_mammals) round(t * (x / sum(x)))
ps_norm_mammals = transform_sample_counts(ps_mammals_sub, standf)

# relative abundances plot - mammals
# NB: to run on server only
load("plot_relab_mammals.RData")

plot_relab_mammals = relab_mammals_plot = ggplot(relab_mammals_df, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("sample") +
  ylab("relative abundance") +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 14)) +
  scale_fill_viridis_d(labels = c("Actinobacteria", "Bacteroidetes", "Cyanobacteria", "Euryarchaeota", "Firmicutes", "Lentisphaerae", "Proteobacteria", "Spirochaetae",  "Tenericutes", "Verrucomicrobia", "Other")) +
  ggtitle("A (mammals)")

plot_relab_mammals

#### _ birds #### 
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
  mutate(RelativeAbundance = Abundance / sum(Abundance))

# relative abundances plot - birds
# can run on local machine
plot_relab_birds = ggplot(relab_birds_df, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("sample") +
  ylab("relative abundance") +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 14)) +
  scale_fill_viridis_d(labels = c("Actinobacteria", "Bacteroidetes", "Cyanobacteria", "Firmicutes", "Lentisphaerae", "Planctomycetes", "Proteobacteria", "Spirochaetae", "Tenericutes", "Verrucomicrobia", "Other")) +
  ggtitle("B (birds)")

plot_relab_birds

#### dominant taxa #### 
# all hosts
fantaxtic::top_taxa(ps_sub, n = 3, tax_level = "phylum")
# Firmicutes 0.395, Proteobacteria 0.316, Bacteroidetes 0.177

fantaxtic::top_taxa(ps_sub, n = 5, tax_level = "class")
# Clostridia 0.233, Gammaproteobacteria 0.169, Bacilli 0.150, Bacteroidia 0.0838, Alphaproteobacteria 0.0714

# mammals
fantaxtic::top_taxa(ps_mammals_sub, n = 3, tax_level = "phylum")
# Firmicutes 0.395, Proteobacteria 0.311, Bacteroidetes 0.184

fantaxtic::top_taxa(ps_mammals_sub, n = 5, tax_level = "class")
# Clostridia 0.245, Gammaproteobacteria 0.168, Bacilli 0.137, Bacteroidia 0.0930, Alphaproteobacteria 0.0711

# birds
fantaxtic::top_taxa(ps_birds_sub, n = 3, tax_level = "phylum")
# Firmicutes 0.396, Proteobacteria 0.338, Bacteroidetes 0.144)

fantaxtic::top_taxa(ps_birds_sub, n = 5, tax_level = "class")
# Bacilli 0.211, Clostridia 0.175, Gammaproteobacteria 0.173, Alphaproteobacteria 0.0730, Betaproteobacteria 0.0645

################ ALPHA DIVERSITY ################

# run PGLMMs on alpha diversity metrics with sociality and diet as fixed effects and host phylogenetic distance as a random effect

#### models ####

# relevel so that "social" is the reference
df_metadata_sub$basic_sociality <- relevel(factor(df_metadata_sub$basic_sociality), ref = "social")

# models for each metric
model_adiv_observed <- pglmm(
  adiv_observed_otus ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

model_adiv_chao1 <- pglmm(
  log(adiv_chao1) ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

model_adiv_shannon <- pglmm(
  log(adiv_shannon) ~ basic_sociality + basic_diet + (1|host_scientific_name) + (1|diet),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

model_adiv_faith <- pglmm(
  log(adiv_faith_pd) ~ basic_sociality + basic_diet + (1|host_scientific_name) + (1|diet),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

# inspect results
summary(model_adiv_observed)
summary(model_adiv_chao1)
summary(model_adiv_shannon)
summary(model_adiv_faith)


# check log-transformed version
model_log <- pglmm(
  log(adiv_observed_otus) ~ basic_sociality + diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
) # results are qualitatively equivalent

#### plots ####

#### _ host class ####
comparisons_class = list(c("c__Aves", "c__Mammalia"))

# obs OTUs
class_obs = df_metadata_sub %>%
  ggplot(aes(x = host_class, y = adiv_observed_otus, fill = host_class)) +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  #xlab("host class") +
  ylab("observed OTUs") +
  scale_x_discrete(labels=c("c__Aves" = "Aves", "c__Mammalia" = "Mammalia")) +
  scale_fill_viridis_d() +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_class, label = "p.signif")

# chao1 
class_chao = df_metadata_sub %>%
  ggplot(aes(x = host_class, y = adiv_chao1, fill = host_class)) +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  #xlab("host class") +
  ylab("Chao1 Index") +
  scale_x_discrete(labels=c("c__Aves" = "Aves", "c__Mammalia" = "Mammalia")) +
  scale_fill_viridis_d() +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_class, label = "p.signif")

# shannon 
class_shannon = df_metadata_sub %>%
  ggplot(aes(x = host_class, y = adiv_shannon, fill = host_class)) +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  xlab("host class") +
  ylab("Shannon Index") +
  scale_x_discrete(labels=c("c__Aves" = "Aves", "c__Mammalia" = "Mammalia")) +
  scale_fill_viridis_d() +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_class, label = "p.signif")

# faith
class_faith = df_metadata_sub %>%
  ggplot(aes(x = host_class, y = adiv_faith_pd, fill = host_class)) +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  xlab("host class") +
  ylab("Faith PD") +
  scale_x_discrete(labels=c("c__Aves" = "Aves", "c__Mammalia" = "Mammalia")) +
  scale_fill_viridis_d() +
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_class, label = "p.signif")

# combine plots in the right order (Observed OTUs, Chao1, Shannon, Faith)
combined_class_plot = (class_obs | class_chao) / (class_shannon | class_faith)
combined_class_plot

#### _ sociality ####
sociality_order = c("solitary", "intermediate", "social")
df_metadata_sub$basic_sociality = factor(df_metadata_sub$basic_sociality, levels = sociality_order)

# obs
soc_obs = df_metadata_sub %>%
  ggplot(aes(x = basic_sociality, y = adiv_observed_otus, fill = basic_sociality)) +
  #geom_violin() +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  xlab("host sociality") +
  ylab("observed OTUs") +
  scale_fill_viridis_d()

# chao1
soc_chao = df_metadata_sub %>%
  ggplot(aes(x = basic_sociality, y = adiv_chao1, fill = basic_sociality)) +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  xlab("host sociality") +
  ylab("Chao1 Index") +
  scale_fill_viridis_d()

# shannon
soc_shannon = df_metadata_sub %>%
  ggplot(aes(x = basic_sociality, y = adiv_shannon, fill = basic_sociality)) +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  xlab("host sociality") +
  ylab("Shannon Index") +
  scale_fill_viridis_d()

# faith
soc_faith = df_metadata_sub %>%
  ggplot(aes(x = basic_sociality, y = adiv_faith_pd, fill = basic_sociality)) +
  #geom_violin() +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  xlab("host sociality") +
  ylab("Faith PD") +
  scale_fill_viridis_d()

# combine sociality plots
combined_soc_plot = (soc_obs | soc_chao) / (soc_shannon | soc_faith)
combined_soc_plot

################ BETA DIVERSITY ################

#### anosim ####

# less informative than db-RDA or PERMANOVA (adonis), because cannot condition on 'random effects', but including it anyway as exploratory
dist_mat <- vegdist(df_otus_sub, method = "robust.aitchison")

# anosim sociality
anosim_result <- anosim(
  dist_mat, 
  df_metadata_sub$basic_sociality, 
  permutations = 999
)

cat("\n=== ANOSIM (alternative test) ===\n")
print(anosim_result)

# anosim diet
anosim_diet <- anosim(
  dist_mat, 
  df_metadata_sub$basic_diet, 
  permutations = 999
)

cat("\n=== ANOSIM (alternative test) ===\n")
print(anosim_diet)

#### db-RDA ####

## model with sociality and diet as fixed effects, study_id and host phylogenetic distance as 'random effects' (in capscale, these are 'Conditions')

# 1) get phylogenetic distances
host_phylo_dist <- cophenetic(host_phylo)

# 2) PCoA to generate phylogenetic axes
# 2a) first, find best k based on scree plot
# create scree plot
host_pcoa <- cmdscale(host_phylo_dist, k = 20, eig = TRUE)  # extract 20 axes initially

# calculate proportion of variance
var_explained <- host_pcoa$eig / sum(host_pcoa$eig)
cumvar <- cumsum(var_explained)

# plot
par(mfrow = c(1, 2))

# scree plot
plot(1:20, var_explained[1:20], 
     type = "b", 
     xlab = "Axis", 
     ylab = "Proportion of Variance",
     main = "Scree Plot",
     pch = 19)
# looks like 5 axes is where it drops off

# cumulative variance
plot(1:20, cumvar[1:20], 
     type = "b",
     xlab = "Axis", 
     ylab = "Cumulative Variance",
     main = "Cumulative Variance Explained",
     pch = 19)
abline(h = 0.80, col = "red", lty = 2)  # 80% threshold; 5 axes explains ~90% variance

par(mfrow = c(1, 1))

# 2b) generate phylogenetic axes
host_pcoa <- cmdscale(host_phylo_dist, k = 5, eig = TRUE)

# 3) match to samples
phylo_coords = as.data.frame(host_pcoa$points)
phylo_coords = phylo_coords %>% 
  rename(PhyPC1 = V1,
         PhyPC2 = V2,
         PhyPC3 = V3,
         PhyPC4 = V4,
         PhyPC5 = V5) %>% 
  rownames_to_column() %>% 
  rename(host_species = rowname)

# 4) add to metadata
df_metadata_sub = df_metadata_sub %>%
  inner_join(phylo_coords, by = "host_species")

# check
df_metadata_sub %>% 
  select(starts_with("PhyPC")) %>% 
  head()

# 5) model with study_id + diet + 5 phylogenetic axes as conditioned terms
mod1_conditioned <- capscale(
  df_otus_sub ~ basic_sociality + 
    Condition(study_id + basic_diet + PhyPC1 + PhyPC2 + PhyPC3 + PhyPC4 + PhyPC5),
  data = df_metadata_sub,
  distance = "robust.aitchison",
  na.action = na.exclude
)

# 6) test model
anova(mod1_conditioned, by = "margin", permutations = 999)

# variance explained
RsquareAdj(mod1_conditioned)

#### _ mammals ####

# add phylogenetic PCoA to metadata
df_metadata_sub_mammals = df_metadata_sub_mammals %>%
  inner_join(phylo_coords, by = "host_species")

# model with study_id + diet + 5 phylogenetic axes as conditioned terms
mod1_conditioned_mammals <- capscale(
  df_otus_sub_mammals ~ basic_sociality + 
    Condition(study_id + basic_diet + PhyPC1 + PhyPC2 + PhyPC3 + PhyPC4 + PhyPC5),
  data = df_metadata_sub_mammals,
  distance = "robust.aitchison",
  na.action = na.exclude
)

# test model
anova(mod1_conditioned_mammals, by = "margin", permutations = 999)

# variance explained
RsquareAdj(mod1_conditioned_mammals)

#### _ birds ####

# add phylogenetic PCoA to metadata
df_metadata_sub_birds = df_metadata_sub_birds %>%
  inner_join(phylo_coords2, by = "host_species")

# model with study_id + diet + 5 phylogenetic axes as conditioned terms
mod1_conditioned_birds <- capscale(
  df_otus_sub_birds ~ basic_sociality + 
    Condition(study_id + basic_diet + PhyPC1 + PhyPC2 + PhyPC3 + PhyPC4 + PhyPC5),
  data = df_metadata_sub_birds,
  distance = "robust.aitchison",
  na.action = na.exclude
)

# test model
anova(mod1_conditioned_birds, by = "margin", permutations = 999)

# variance explained
RsquareAdj(mod1_conditioned_birds)


#### adonis ####
mod1_adonis <- adonis2(
  dist_mat ~ basic_sociality + study_id + basic_diet + 
    PhyPC1 + PhyPC2 + PhyPC3 + PhyPC4 + PhyPC5,
  data = df_metadata_sub,
  permutations = 999,
  by = "margin"  # tests each term after all others
)

mod1_adonis

#### betadisp ####

# measures dispersion within each group and compared to dispersion between group
# calculate distances
dist_mat <- vegdist(df_otus_sub, method = "robust.aitchison")

# test dispersion differences among sociality groups
disp_test_social <- betadisper(dist_mat, df_metadata_sub$basic_sociality)
anova(disp_test_social) # p < 0.001

plot(disp_test_social)
boxplot(disp_test_social, main = "Distance to centroid by social behaviour")

cat("\nPairwise dispersion differences:\n")
print(TukeyHSD(disp_test_social))

# test dietary groups
disp_test_diet <- betadisper(dist_mat, df_metadata_sub$basic_diet)
anova(disp_test_diet) # p < 0.001

plot(disp_test_diet)
boxplot(disp_test_diet, main = "Distance to centroid by diet")
# calculate distances
dist_mat <- vegdist(df_otus_sub, method = "robust.aitchison")

# test dispersion differences among sociality groups
disp_test_social <- betadisper(dist_mat, df_metadata_sub$basic_sociality)
anova(disp_test_social) # p < 0.001

plot(disp_test_social)
boxplot(disp_test_social, main = "Distance to centroid by social behaviour")

cat("\nPairwise dispersion differences:\n")
print(TukeyHSD(disp_test_social))

# test dietary groups
disp_test_diet <- betadisper(dist_mat, df_metadata_sub$basic_diet)
anova(disp_test_diet) # p < 0.001

plot(disp_test_diet)
boxplot(disp_test_diet, main = "Distance to centroid by diet")


# _ plot all hosts ----
df_metadata_sub = df_metadata_sub %>% rownames_to_column() %>% rename("sample_id" = "rowname")

# species
rda_scores_sp_df = rda_sp_all %>% 
  scores(display = "sites") %>% # extract the site scores
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(sample_id = rowname) %>% #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_sp_all, groups = df_metadata_sub$host_species, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 4.35%, CAP2 1.39%

rda_scores_sp_df %>% 
  ggplot(aes(x = CAP1, y = CAP2, colour = host_species, shape = host_class)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  scale_shape_manual(values = c("c__Aves" = 16, "c__Mammalia" = 17), labels = c("Aves", "Mammalia")) +
  labs(
    colour = "host species", 
    shape = "host class",
    x = "CAP1 (4.35%)", 
    y = "CAP2 (1.39%)"
  ) +
  theme(text = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  guides(colour = "none") +
  ggtitle("A")

# diet
rda_scores_diet_df = rda_diet_all %>% 
  scores(display = "sites") %>% # extract the site scores
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(sample_id = rowname) %>% #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_diet_all, groups = df_metadata_sub$basic_diet, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 1.13%, CAP2 0.67%

rda_scores_diet_df %>% 
  ggplot(aes(x = CAP1, y = CAP2, colour = basic_diet, shape = host_class)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  scale_shape_manual(values = c("c__Aves" = 16, "c__Mammalia" = 17), labels = c("Aves", "Mammalia")) +
  labs(
    colour = "diet", 
    shape = "host class",
    x = "CAP1 (1.13%)", 
    y = "CAP2 (0.67%)"
  ) +
  theme(text = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  ggtitle("B")

# sociality
rda_scores_social_df = rda_social_all %>% 
  scores(display = "sites") %>% # extract the site scores
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(sample_id = rowname) %>% #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_social_all, groups = df_metadata_sub$basic_sociality, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 1.66%, CAP2 0.52%

rda_scores_social_df %>% 
  ggplot(aes(x = CAP1, y = CAP2, colour = basic_sociality, shape = host_class)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  scale_shape_manual(values = c("c__Aves" = 16, "c__Mammalia" = 17), labels = c("Aves", "Mammalia")) +
  labs(
    colour = "sociality", 
    shape = "order",
    x = "CAP1 (1.66%)", 
    y = "CAP2 (0.52%)"
  ) +
  theme(text = element_text(size = 14),
        legend.text = element_text(size = 14)) +
  ggtitle("C")

# _ plot mammals ----
# species
rda_scores_sp_mammals_df = rda_sp_mammals %>% 
  scores(display = "sites") %>% # extract the site scores
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(sample_id = rowname) %>% #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_sp_mammals, groups = df_metadata_sub_mammals$host_species, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 4.52%, CAP2 1.57%

rda_scores_sp_mammals_df %>% 
  ggplot(aes(x = CAP1, y = CAP2, colour = host_species, shape = basic_diet)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  labs(
    colour = "host species", 
    shape = "diet",
    x = "CAP1 (4.52%)", 
    y = "CAP2 (1.57%)"
  ) +
  theme(
    text = element_text(size = 14), 
    legend.position = "none") +
  ggtitle("A (mammals)")

# diet
rda_scores_diet_mammals_df = rda_diet_mammals %>% 
  scores(display = "sites") %>% # extract the site scores
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(sample_id = rowname) %>% #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_diet_mammals, groups = df_metadata_sub_mammals$basic_diet, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 1.28%, CAP2 0.19%

rda_scores_diet_mammals_df %>% 
  ggplot(aes(x = CAP1, y = CAP2, colour = basic_diet, shape = basic_sociality)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  labs(
    colour = "diet", 
    shape = "sociality",
    x = "CAP1 (1.28%)", 
    y = "CAP2 (0.19%)"
  ) +
  theme(
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  ggtitle("B (mammals)")

# sociality
rda_scores_social_mammals_df = rda_social_mammals %>% 
  scores(display = "sites") %>% # extract the site scores
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(sample_id = rowname) %>% #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_social_mammals, groups = df_metadata_sub_mammals$basic_sociality, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 2.03%, CAP2 0.70%

rda_scores_social_mammals_df %>% 
  ggplot(aes(x = CAP1, y = CAP2, colour = basic_sociality)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  labs(
    colour = "sociality", 
    x = "CAP1 (2.03%)", 
    y = "CAP2 (0.70%)"
  ) +
  theme(
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  ggtitle("C (mammals)")

# _ plot birds ----
df_metadata_sub_birds = df_metadata_sub_birds %>% rownames_to_column() %>% rename("sample_id" = "rowname")

# species
rda_scores_sp_birds_df = rda_sp_birds %>% 
  scores(display = "sites") %>% # extract the site scores
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(sample_id = rowname) %>% #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_sp_birds, groups = df_metadata_sub_birds$host_species, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 10.65%, CAP2 4.17%

rda_scores_sp_birds_df %>% 
  ggplot(aes(x = CAP1, y = CAP2, colour = host_scientific_name, shape = basic_diet)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  labs(
    colour = "host species", 
    shape = "diet",
    x = "CAP1 (10.65%)", 
    y = "CAP2 (4.17%)"
  ) +
  theme(
    text = element_text(size = 14),
    legend.text = element_text(face = "italic")) +
  ggtitle("A (birds)")

# diet
rda_scores_diet_birds_df = rda_diet_birds %>% 
  scores(display = "sites") %>% # extract the site scores
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(sample_id = rowname) %>% # add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_diet_birds, groups = df_metadata_sub_birds$basic_diet, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 3.82% (no other axes)

rda_scores_diet_birds_df %>% 
  ggplot(aes(x = CAP1, y = MDS1, colour = basic_diet, shape = basic_sociality)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  labs(
    colour = "diet", 
    shape = "sociality",
    x = "CAP1 (1.28%)", 
    y = "CAP2 (0.19%)"
  ) +
  theme(
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  ggtitle("B (birds)")

# sociality
rda_scores_social_birds_df = rda_social_birds %>% 
  scores(display = "sites") %>% # extract the site scores
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(sample_id = rowname) %>% #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_social_birds, groups = df_metadata_sub_birds$basic_sociality, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 4.07%, CAP2 1.16%

rda_scores_social_birds_df %>% 
  ggplot(aes(x = CAP1, y = CAP2, colour = basic_sociality)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dotted") + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_colour_viridis_d() +
  labs(
    colour = "sociality", 
    x = "CAP1 (4.07%)", 
    y = "CAP2 (1.16%)"
  ) +
  theme(
    legend.text = element_text(size = 14), 
    legend.title = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  ggtitle("C (birds)")

################ BETADISP ################

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

# check model assumptions
plot(lm_weighted) # residuals vs fitted (homoscedasticity) and qq plot (normality of residuals) look good

#### plot #####

# ordering
disp_data = disp_data %>% 
  mutate(sociality = fct_relevel(sociality, c("solitary", "intermediate", "social"))) %>%
  arrange(sociality, n) %>%  # reorder by sociality then sample size
  mutate(host_scientific_name = factor(host_scientific_name, levels = unique(host_scientific_name)))

# plot dispersion distances
disp_data %>%
  ggplot(aes(x = host_scientific_name, y = distance, fill = sociality)) +
  geom_boxplot(width = 0.6) +
  scale_fill_viridis_d(name = "Sociality") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 140, by = 20)) +
  geom_text(data = disp_data, aes(host_scientific_name, Inf, label = n), hjust = "inward") +
  labs(x = "Host species", y = "Distance from centroid") +
  theme(text = element_text(size = 14),
        axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip()

################ DIFFERENTIAL ABUNDANCE ################

##### setup #####

# all hosts
df_metadata_sub_da = df_metadata_sub %>% 
  filter(basic_sociality != "intermediate") %>% # exclude intermediate sociality
  column_to_rownames("sample_id") # bring rownames back

# subsample otu table
common_samples_sub_da = intersect(rownames(mx_otus), rownames(df_metadata_sub_da))

df_otus_sub_da = df_otus_sub %>%
  rownames_to_column() %>%
  filter(rowname %in% common_samples_sub_da) %>% 
  column_to_rownames("rowname") %>%
  select(where(~ sum(.) != 0))

# mammals
df_metadata_sub_mammals_da = df_metadata_sub %>% 
  filter(basic_sociality != "intermediate") %>%
  filter(host_class == "c__Mammalia") %>% 
  column_to_rownames("sample_id")

common_samples_sub_mammals_da = intersect(rownames(mx_otus), rownames(df_metadata_sub_mammals_da))

df_otus_sub_mammals_da = df_otus_sub %>%
  rownames_to_column() %>%
  filter(rowname %in% common_samples_sub_mammals_da) %>% 
  column_to_rownames("rowname") %>%
  select(where(~ sum(.) != 0))

# birds
df_metadata_sub_birds_da = df_metadata_sub %>% 
  filter(basic_sociality != "intermediate") %>%
  filter(host_class == "c__Aves") %>% 
  column_to_rownames("sample_id")

common_samples_sub_birds_da = intersect(rownames(mx_otus), rownames(df_metadata_sub_birds_da))

df_otus_sub_birds_da = df_otus_sub %>%
  rownames_to_column() %>%
  filter(rowname %in% common_samples_sub_birds_da) %>% 
  column_to_rownames("rowname") %>%
  select(where(~ sum(.) != 0))

# run this on the server only
# all hosts
v_social_da = df_metadata_sub_da %>% pull(basic_sociality)
social_clr = aldex.clr(t(df_otus_sub_da), v_social_da, mc.samples = 200) #basic_sociality is a character vector, so aldex sorts it alphabetically -> social is the reference group
social_ttest = aldex.ttest(social_clr) #do not run on laptop
aldex_social_effect = aldex.effect(social_clr, CI = TRUE)
social_aldex_all = data.frame(social_ttest, aldex_social_effect)
social_aldex_all = social_aldex_all %>% rownames_to_column()

# mammals
v_social_da_mammals = df_metadata_sub_mammals_da %>% pull(basic_sociality)
social_clr_mammals = aldex.clr(t(df_otus_sub_mammals_da), v_social_da_mammals, mc.samples = 200)
social_ttest_mammals = aldex.ttest(social_clr_mammals)
aldex_social_effect_mammals = aldex.effect(social_clr_mammals, CI = TRUE)
social_aldex_mammals = data.frame(social_ttest_mammals, aldex_social_effect_mammals)
social_aldex_mammals = social_aldex_mammals %>% rownames_to_column()

# birds
v_social_da_birds = df_metadata_sub_birds_da %>% pull(basic_sociality)
social_clr_birds = aldex.clr(t(df_otus_sub_birds_da), v_social_da_birds, mc.samples = 200)
social_ttest_birds = aldex.ttest(social_clr_birds)
aldex_social_effect_birds = aldex.effect(social_clr_birds, CI = TRUE)
social_aldex_birds = data.frame(social_ttest_birds, aldex_social_effect_birds)
social_aldex_birds = social_aldex_birds %>% rownames_to_column()

write_csv(social_aldex_all, "social_aldex_all.csv")
write_csv(social_aldex_mammals, "social_aldex_mammals.csv")
write_csv(social_aldex_birds, "social_aldex_birds.csv")

save(social_aldex_all, social_aldex_mammals, social_aldex_birds, file = "diff_ab_data.RData")

#### data input #####
social_aldex_all = read_csv("social_aldex_all.csv")
social_aldex_mammals = read_csv("social_aldex_mammals.csv")
social_aldex_birds = read_csv("social_aldex_birds.csv")

#### plots #####
par(mfrow = c(1,2))
aldex.plot(social_aldex_all, type = "MA", test = "welch", main = "MA plot")
aldex.plot(social_aldex_all, type = "MW", test = "welch", main = "effect plot")

# identify OTUs by name
taxonomy_table = p_tax %>% as.data.frame() %>% rownames_to_column()

# generate df of differentially abundant OTUs
social_diff_asvs_mammals = social_aldex_mammals %>% 
  filter(we.eBH < 0.001) %>% 
  rownames_to_column("...1") %>% 
  left_join(taxonomy_table, by = "rowname") %>% # add ASV taxonomical info
  replace_na(list(genus = "")) %>%
  replace_na(list(species = "sp.")) %>%
  mutate(across(c("genus"), substr, 6, nchar(genus))) %>% 
  mutate(across(c("species"), substr, 6, nchar(species))) %>% 
  mutate(ci_product = effect.low * effect.high ) %>% #column that returns TRUE if CI does not contain 0# %>% 
  mutate(ci_no_zero = ci_product > 0) #returns TRUE if CI does not include zero

social_diff_asvs_mammals %>% 
  ggplot(aes(x = effect, y = reorder(species, (effect*-1)), colour = phylum)) +
  geom_point(size = 2.5) +
  scale_colour_viridis_d(labels = c("Actinobacteria", "Bacteroidetes", "Firmicutes", "Proteobacteria")) +
  xlab("effect size") +
  ylab("OTU") +
  theme(axis.text.y = element_text(face = "italic"),
        text = element_text(size = 14)) +
  #ggtitle("Differentially abundant OTUs in mammals,\nsocial v solitary\n(BH-corrected)") + 
  geom_vline(xintercept = 0, linetype = "dotted")
# since social is the reference group, effect < 0 means OTU is higher in social

# birds
social_diff_asvs_birds = social_aldex_birds %>%
  filter(we.eBH < 0.05) %>% #we.eBH = expected Benjamini-Hochberg corrected p-value of Welch’s t test (higher p still yields 0 differentially abundant OTUs)
  rownames_to_column("...1") %>% 
  left_join(taxonomy_table, by = "rowname") %>% # add ASV taxonomical info
  replace_na(list(genus = "")) %>%
  replace_na(list(species = "sp.")) %>%
  mutate(otu_scientific_name = paste(genus, species, sep = " ")) %>% #concatenate to full ASV name
  mutate(otu_scientific_name = str_trim(otu_scientific_name, side = "left")) %>%
  mutate(otu_scientific_name = as.factor(otu_scientific_name)) %>% 
  mutate(ci_product = effect.low * effect.high ) %>% #column that returns TRUE if CI does not contain 0# %>% 
  mutate(ci_no_zero = ci_product > 0) #returns TRUE if CI does not include zero
#there are no differentially abundant OTUs in birds