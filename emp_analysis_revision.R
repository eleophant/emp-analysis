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

load(file = "emp_analysis_data_2025_02_03.RData") #this dataset has been modified to omit 1 reptile species, black bear, colobine primates

################ SUBSAMPLING ################

#### make dfs ####
df_obs_n = df_metadata |> 
  group_by(host_species) |> 
  summarise(n = n()) |> 
  print(n = 27)

df_obs_n |> summarise(mean(n), median(n), max(n), min(n))
#mean 41, median 3, max 902, min 1

# subsample metadata
df_metadata_sub = df_metadata |> 
  rownames_to_column() |> 
  group_by(host_species) |> 
  slice_sample(n = 41)  |>
  ungroup() |> 
  column_to_rownames("rowname")

# subsample otu table
common_samples_sub = intersect(rownames(mx_otus), rownames(df_metadata_sub))

df_otus_sub = mx_otus |> as.data.frame() |> 
  rownames_to_column() |> 
  filter(rowname %in% common_samples_sub) |> 
  column_to_rownames("rowname") |>
  select(where(~ sum(.) != 0)) #remove empty columns (OTUs that are 0 everywhere)

# ensure that row order matches
df_metadata_sub = df_metadata_sub |>
  slice(match(rownames(df_otus_sub), rownames(df_metadata_sub)))
all(rownames(df_otus_sub) == rownames(df_metadata_sub)) #TRUE

#### phyloseq ####
# set up metadata & otu table
p_metadata_sub = sample_data(df_metadata_sub)
p_otu_sub = otu_table(df_otus_sub, taxa_are_rows = FALSE)

# check that sample names match
all(sample_names(p_metadata_sub) == sample_names(p_otu_sub)) #TRUE

# assemble phyloseq object
physeq_sub = phyloseq(p_otu_sub, p_metadata_sub, p_tax)
ps_sub = merge_phyloseq(physeq_sub, tree_emp)

#### subset classes ####
# mammals
df_metadata_sub_mammals = df_metadata_sub |>
  filter(host_class == "c__Mammalia")
common_samples_sub_mammals = intersect(rownames(mx_otus), rownames(df_metadata_sub_mammals))

df_otus_sub_mammals = mx_otus |> as.data.frame() |> 
  rownames_to_column() |> 
  filter(rowname %in% common_samples_sub_mammals) |> 
  column_to_rownames("rowname") |>
  select(where(~ sum(.) != 0))

# ensure that row order matches
all(rownames(df_otus_sub_mammals) == rownames(df_metadata_sub_mammals)) #TRUE

# birds
df_metadata_sub_birds = df_metadata_sub |>
  filter(host_class == "c__Aves")
common_samples_sub_birds = intersect(rownames(mx_otus), rownames(df_metadata_sub_birds))

df_otus_sub_birds = mx_otus |> as.data.frame() |> 
  rownames_to_column() |> 
  filter(rowname %in% common_samples_sub_birds) |> 
  column_to_rownames("rowname") |>
  select(where(~ sum(.) != 0))

# ensure that row order matches
df_metadata_sub_birds = df_metadata_sub_birds |> 
  slice(match(rownames(df_otus_sub_birds), rownames(df_metadata_sub_birds)))
all(rownames(df_otus_sub_birds) == rownames(df_metadata_sub_birds)) #TRUE

# phyloseq
ps_mammals_sub = subset_samples(ps_sub, host_class == "c__Mammalia")
ps_birds_sub = subset_samples(ps_sub, host_class == "c__Aves")

################ DATA EXPLORATION ################

#### all hosts ####
df_metadata_sub |> nrow() #375 samples
df_metadata_sub |> group_by(host_species) |> 
  summarise(n = n()) |> print(n = 100) #across 45 species
df_otus_sub |> as.data.frame() |> ncol() #number of OTUs = 34,482
df_otus_sub |> sum() #total reads = 13,085,267 (this number will vary slightly each time this script is run since subsampling is random)

# per sample
mean(rowSums(df_otus_sub)) #mean number of reads = 34,894.05
sd(rowSums(df_otus_sub)) #24,286.86

otus_per_sample = df_otus_sub |> mutate(total_otus = rowSums(df_otus_sub != 0)) |> rownames_to_column() |> 
  select(rowname, total_otus)

otus_per_sample = df_metadata_sub |> 
  rownames_to_column() |> 
  #rename(rowname = sample_id) |> 
  full_join(otus_per_sample) |> 
  select(rowname, host_class, host_scientific_name, total_otus)

otus_per_sample |> summarise(mean(total_otus), median(total_otus), max(total_otus), min(total_otus))
#mean 1019, median 961, max 2817, min 114

#### mammals ####
# collectively
df_metadata_sub_mammals |> nrow() #312 samples
df_metadata_sub_mammals |> group_by(host_species) |> 
  summarise(n = n()) |> print(n = 35) #across 35 species
df_otus_sub_mammals |> ncol() #31,451 OTUs
df_otus_sub_mammals |> sum() #total reads = 11,086,194 (this number will vary slightly each time this script is run since subsampling is random)

# per sample
mean(rowSums(df_otus_sub_mammals)) #mean number of reads = 35,532.67
sd(rowSums(df_otus_sub_mammals)) #23,866.22

otus_per_sample_mammals = df_otus_sub_mammals |> mutate(total_otus = rowSums(df_otus_sub_mammals != 0)) |> rownames_to_column() |> 
  select(rowname, total_otus)

otus_per_sample_mammals = df_metadata_sub_mammals |> 
  #rename(rowname = sample_id) |> 
  rownames_to_column() |> 
  full_join(otus_per_sample_mammals) |> 
  select(rowname, host_class, host_scientific_name, total_otus)

otus_per_sample_mammals |> summarise(mean(total_otus), median(total_otus), max(total_otus), min(total_otus))
#mean 1062, median 1004, max 2473, min 209

#### birds ####
df_metadata_sub_birds |> nrow() #63 samples
df_metadata_sub_birds |> group_by(host_species) |> 
  summarise(n = n()) |> print(n = 100) #across 10 species
df_otus_sub_birds |> ncol() #13,915
df_otus_sub_birds |> sum() #total reads = 1,999,073 (this number will vary slightly each time this script is run since subsampling is random)

# per sample
mean(rowSums(df_otus_sub_birds)) #mean number of reads = 31,731.32
sd(rowSums(df_otus_sub_birds)) #26,244.05

# per sample
otus_per_sample_birds = df_otus_sub_birds |> mutate(total_otus = rowSums(df_otus_sub_birds != 0)) |> rownames_to_column() |> 
  select(rowname, total_otus)

otus_per_sample_birds = df_metadata_sub_birds |> 
  #rename(rowname = sample_id) |> 
  rownames_to_column() |> 
  full_join(otus_per_sample_birds) |> 
  select(rowname, host_class, host_scientific_name, total_otus)

otus_per_sample_birds |> summarise(mean(total_otus), median(total_otus), max(total_otus), min(total_otus))
#mean 807, median 695, max 2817, min 114

#### compare classes ####
wilcox.test(total_otus ~ host_class, data = otus_per_sample) #W = 6841, p-value = 0.0001415

#### social OTUs ####
# goal: find (Limosi)lactobacillus and Bifidobacterium in df_otus
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
# Bifido: 2685/13085267*100 = 0.02051926 %

# Limosilact OTUs
df_otus_sub %>% 
  select(all_of(ls_limosi)) %>% # all Limosilact OTUs
  summarise(across(everything(), ~ sum(., na.rm = TRUE))) %>% 
  sum() # 59335
# Limosilact: 59335/13085267*100 = 0.4534489 %

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
relab_birds_df = ps_norm_birds |>
  psmelt() |>
  group_by(Sample, phylum) |>
  summarize(Abundance = sum(Abundance), .groups = "drop") |>
  mutate(RelativeAbundance = Abundance / sum(Abundance))

# get top 10 phyla
top_phyla_birds = relab_birds_df |>
  group_by(phylum) |>
  summarize(TotalAbundance = sum(Abundance), .groups = "drop") |>
  arrange(desc(TotalAbundance)) |>
  slice_head(n = 10) |>
  pull(phylum)

# reclassify other phyla as "Other"
relab_birds_df = relab_birds_df |>
  mutate(phylum = ifelse(phylum %in% top_phyla_birds, phylum, "Other"))

# recalculate relative abundances with new categories
relab_birds_df = relab_birds_df |>
  group_by(Sample, phylum) |>
  summarize(Abundance = sum(Abundance), .groups = "drop") |>
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
#Firmicutes 0.395, Proteobacteria 0.316, Bacteroidetes 0.177

fantaxtic::top_taxa(ps_sub, n = 5, tax_level = "class")
#Clostridia 0.233, Gammaproteobacteria 0.169, Bacilli 0.150, Bacteroidia 0.0838, Alphaproteobacteria 0.0714

# mammals
fantaxtic::top_taxa(ps_mammals_sub, n = 3, tax_level = "phylum")
#Firmicutes 0.395, Proteobacteria 0.311, Bacteroidetes 0.184

fantaxtic::top_taxa(ps_mammals_sub, n = 5, tax_level = "class")
#Clostridia 0.245, Gammaproteobacteria 0.168, Bacilli 0.137, Bacteroidia 0.0930, Alphaproteobacteria 0.0711

# birds
fantaxtic::top_taxa(ps_birds_sub, n = 3, tax_level = "phylum")
#Firmicutes 0.396, Proteobacteria 0.338, Bacteroidetes 0.144)

fantaxtic::top_taxa(ps_birds_sub, n = 5, tax_level = "class")
#Bacilli 0.211, Clostridia 0.175, Gammaproteobacteria 0.173, Alphaproteobacteria 0.0730, Betaproteobacteria 0.0645

################ HOST PHYLOGENY ################


# fix metadata to match OTL
df_metadata_tree <- df_metadata_sub %>%
  # remove "Genus sp."
  filter(host_scientific_name != "Anser sp." & 
           host_scientific_name != "Macropus sp.") %>% 
  # fix synonyms 
  mutate(
    host_scientific_name = case_when(
      host_scientific_name == "Proteles cristatus" ~ "Proteles cristata",
      host_scientific_name == "Thraupis palmarum" ~ "Tangara palmarum",
      TRUE ~ host_scientific_name
    )
  )

# get OTT IDs for all EMP hosts
ott_ids <- df_metadata_tree %>%
  pull(host_scientific_name) %>% # extract host names
  unique() %>% # remove duplicates
  tnrs_match_names() %>% # match to OpenTree
  pull(ott_id) # extract OTT IDs

# build tree from species in df_metadata_tree
host_phylo <- tol_induced_subtree(ott_ids = ott_ids) %>%
  compute.brlen()

# plot & check tip labels
ape::plot.phylo(host_phylo, cex = 0.6)
head(host_phylo$tip.label)

# remove OTT suffix from tip labels
host_phylo$tip.label = host_phylo$tip.label %>%
  gsub("_ott[0-9]+$", "", .) %>% # remove "_ottXXXX" suffix
  gsub("_", " ", .) %>% # replace underscores with spaces
  trimws() # remove any stray whitespaces 

# re-check tip labels
head(host_phylo$tip.label) # now just "Genus species"
ape::plot.phylo(host_phylo, cex = 0.6)

all(df_metadata_tree$host_scientific_name %in% host_phylo$tip.label) #TRUE

################ ALPHA DIVERSITY ################

# run PGLMMs on alpha diversity metrics with sociality and diet as fixed effects and host phylogenetic distance as a random effect

#### models ####

# models for each metric
model_adiv_observed <- pglmm(
  adiv_observed_otus ~ basic_sociality + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

model_adiv_chao1 <- pglmm(
  adiv_chao1 ~ basic_sociality + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

model_adiv_shannon <- pglmm(
  adiv_shannon ~ basic_sociality + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

model_adiv_faith <- pglmm(
  adiv_faith_pd ~ basic_sociality + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

# inspect results
summary(model_adiv_observed)
summary(model_adiv_chao1)
summary(model_adiv_shannon)
summary(model_adiv_faith)

#### plots ####

#### _ host class ####
comparisons_class = list(c("c__Aves", "c__Mammalia"))

# obs OTUs
class_obs = df_metadata_sub |>
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
class_chao = df_metadata_sub |>
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
class_shannon = df_metadata_sub |>
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
class_faith = df_metadata_sub |>
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
soc_obs = df_metadata_sub |>
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
soc_chao = df_metadata_sub |>
  ggplot(aes(x = basic_sociality, y = adiv_chao1, fill = basic_sociality)) +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(), axis.text.x = element_blank()) +
  xlab("host sociality") +
  ylab("Chao1 Index") +
  scale_fill_viridis_d()

# shannon
soc_shannon = df_metadata_sub |>
  ggplot(aes(x = basic_sociality, y = adiv_shannon, fill = basic_sociality)) +
  geom_boxplot() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14), axis.text = element_text(size = 14)) +
  xlab("host sociality") +
  ylab("Shannon Index") +
  scale_fill_viridis_d()

# faith
soc_faith = df_metadata_sub |>
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

################ DISSIMILARITY ################

#### mammals ####

aitchison_mammals = vegdist(df_otus_sub_mammals[,-1], method = "robust.aitchison", pseudocount = 1) #aitchison performs a clr, which only works on positive values: need to assign 1 to counts of 0

# make a table with 3 columns: column id, row id, dissimilarity value
dissimilarity_table_mammals = bind_cols(t(combn(nrow(df_otus_sub_mammals), 2)), diss = aitchison_mammals)

# assign row numbers to join to dissimilarity table
df_metadata_sub_mammals = df_metadata_sub_mammals |> mutate(row_number = row_number())

dissimilarity_table_mammals = dissimilarity_table_mammals |>
  rename(row_number = ...1, column_number = ...2) |>
  left_join(select(df_metadata_sub_mammals, row_number, host_scientific_name, basic_sociality), by = "row_number") |>
  rename(host_species_1 = host_scientific_name)

dissimilarity_table_mammals_sp2 = dissimilarity_table_mammals |> 
  select(column_number, diss) |> 
  left_join(select(df_metadata_sub_mammals, host_scientific_name, row_number), by = c("column_number" = "row_number")) |> rename(host_species_2 = host_scientific_name)

dissimilarity_table_mammals = dissimilarity_table_mammals |> 
  full_join(dissimilarity_table_mammals_sp2)

#### birds ####

aitchison_birds = vegdist(df_otus_sub_birds[,-1], method = "robust.aitchison", pseudocount = 1)

# make a table with 3 columns: column id, row id, dissimilarity value
dissimilarity_table_birds = bind_cols(t(combn(nrow(df_otus_sub_birds), 2)), diss = aitchison_birds)

# assign row numbers to join to dissimilarity table
df_metadata_sub_birds = df_metadata_sub_birds |> mutate(row_number = row_number())

dissimilarity_table_birds = dissimilarity_table_birds |>
  rename(row_number = ...1, column_number = ...2) |>
  left_join(select(df_metadata_sub_birds, row_number, host_scientific_name, basic_sociality), by = "row_number") |>
  rename(host_species_1 = host_scientific_name)

dissimilarity_table_birds_sp2 = dissimilarity_table_birds |> 
  select(column_number, diss) |> 
  left_join(select(df_metadata_sub_birds, host_scientific_name, row_number), by = c("column_number" = "row_number")) |> rename(host_species_2 = host_scientific_name)

dissimilarity_table_birds = dissimilarity_table_birds |> 
  full_join(dissimilarity_table_birds_sp2)

#### plots ####

# mammals
# calculate sample sizes per species
sample_size_mammals = df_metadata_sub_mammals  |> 
  group_by(host_scientific_name, basic_sociality) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(basic_sociality = fct_relevel(basic_sociality, c("social", "intermediate", "solitary"))) |>
  arrange(basic_sociality, n) |>  # reorder species first by sociality, then by sample size within each group
  mutate(host_scientific_name = factor(host_scientific_name, levels = unique(host_scientific_name))) |>
  rename(host_species_1 = host_scientific_name) |>
  filter(n != 1) 

diss_mammals = sample_size_mammals |> pull(host_species_1)

# plot dissimilarity values
dissimilarity_table_mammals |>
  filter(host_species_1 %in% diss_mammals) |>
  filter(host_species_1 == host_species_2) |>
  mutate(
    basic_sociality = fct_relevel(basic_sociality, "social", "intermediate", "solitary"),
    host_species_1 = factor(host_species_1, levels = levels(sample_size_mammals$host_species_1)) # ensure order is consistent
  ) |>
  ggplot(aes(x = host_species_1, y = diss, fill = basic_sociality)) +
  geom_boxplot(width = 0.8) +
  scale_color_viridis(discrete = TRUE, aesthetics = "fill", name = "sociality") +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0, 140, by = 40)) +
  geom_text(data = sample_size_mammals, aes(host_species_1, Inf, label = n), hjust = "inward") +
  labs(x = "host species", y = "dissimilarity (robust Aitchison)") +
  theme(text = element_text(size = 14),
        axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip()

# birds 
# calculate sample sizes per species
sample_size_birds = df_metadata_sub_birds  |> 
  group_by(host_scientific_name, basic_sociality) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(basic_sociality = fct_relevel(basic_sociality, c("social", "intermediate", "solitary"))) |>
  arrange(basic_sociality, n) |>  # reorder species first by sociality, then by sample size within each group
  mutate(host_scientific_name = factor(host_scientific_name, levels = unique(host_scientific_name))) |>
  rename(host_species_1 = host_scientific_name) |>
  filter(n != 1) 

diss_birds = sample_size_birds |> pull(host_species_1)

# plot dissimilarity values
dissimilarity_table_birds |>
  filter(host_species_1 %in% diss_birds) |>
  filter(host_species_1 == host_species_2) |>
  mutate(
    basic_sociality = fct_relevel(basic_sociality, "social", "intermediate", "solitary"),
    host_species_1 = factor(host_species_1, levels = levels(sample_size_birds$host_species_1)) # ensure order is consistent
  ) |>
  ggplot(aes(x = host_species_1, y = diss, fill = basic_sociality)) +
  geom_boxplot(width = 0.8) +
  scale_color_viridis(discrete = TRUE, aesthetics = "fill", name = "sociality") +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0, 140, by = 40)) +
  geom_text(data = sample_size_birds, aes(host_species_1, Inf, label = n), hjust = "inward") +
  labs(x = "host species", y = "dissimilarity (robust Aitchison)") +
  theme(text = element_text(size = 14),
        axis.text.y = element_text(face = "italic")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  coord_flip()

#### stats ####
# run these on server only
glm_diss_mammals = lmer(diss ~ basic_sociality + (1|host_species_1), data = dissimilarity_table_mammals)
summary(glm_diss_mammals)
jtools::summ(glm_diss_mammals) #Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()
#server results solitary β = -0.41, SE = 3.91, t = -0.11, p = 0.39; social: β = -3.25, SE = 3.72, t = -0.87, p = 0.92

glm_diss_birds = lmer(diss ~ basic_sociality + (1|host_species_1), data = dissimilarity_table_birds)
summary(glm_diss_birds)
jtools::summ(glm_diss_birds) #social p = 0.60, solitary 0.46

################ DB-RDA ################

#### all hosts ####

# ensure that row order matches
all(rownames(df_otus_sub) == rownames(df_metadata_sub)) #TRUE

# by species
rda_sp_all = capscale(formula = df_otus_sub ~ host_species, data = df_metadata_sub,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_sp_all = anova(rda_sp_all) 
print(aov_rda_sp_all) #p = 0.001
RsquareAdj(rda_sp_all) #adj R^2 = 0.08435098

# by diet
rda_diet_all = capscale(formula = df_otus_sub ~ basic_diet, data = df_metadata_sub,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_diet_all = anova(rda_diet_all)
print(aov_rda_diet_all) #p = 0.001
RsquareAdj(rda_diet_all) #adj R^2 = 0.01276794

# by sociality
rda_social_all = capscale(formula = df_otus_sub ~ basic_sociality, data = df_metadata_sub,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_social_all = anova(rda_social_all)
print(aov_rda_diet_all) # p = 0.001
RsquareAdj(rda_social_all) #ajd R^2 = 0.01653155

#### full model ####
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
host_phylo_dist <- cophenetic(host_tree)
host_pcoa <- cmdscale(host_phylo_dist, k = 5, eig = TRUE)

# 3) match to samples
phylo_coords <- as.data.frame(host_pcoa$points)
colnames(phylo_coords) <- paste0("PhyPC", 1:5)
species_vector <- as.character(df_metadata_sub$host_scientific_name)
sample_phylo <- phylo_coords[species_vector, ]
rownames(sample_phylo) <- rownames(df_metadata_sub)

# 4) add to metadata
df_metadata_sub <- cbind(df_metadata_sub, sample_phylo)

# 5) model with study_id and 5 phylogenetic axes as conditioning terms
mod1_conditioned <- capscale(
  formula = df_otus_sub ~ basic_diet + basic_sociality + 
    Condition(study_id + PhyPC1 + PhyPC2 + PhyPC3 + PhyPC4 + PhyPC5),
  data = df_metadata_sub,
  distance = "robust.aitchison",
  na.action = na.exclude
)

# 6) test model
anova(mod1_conditioned, by = "margin", permutations = 999)

# 6. Variance partitioning
summary(mod1_conditioned)

############ PREVIOUS CODE #############
# run null and full models
mod0_all = capscale(df_otus_sub ~ 1, data = df_metadata_sub, distance = "robust.aitchison", na.action = na.exclude)
mod1_all = capscale(formula = df_otus_sub ~ host_species + basic_diet + basic_sociality + (1|study_id), data = df_metadata_sub,  distance = "robust.aitchison", na.action = na.exclude) #Some constraints or conditions were aliased because they were redundant. This can happen if terms are linearly dependent (collinear): ‘basic_dietherbivorous’, 'basic_dietomnivorous’, ‘basic_socialitysolitary’

# ordistep
step_all = ordistep(mod0_all, scope = formula(mod1_all), direction = "both", perm.max = 200, na.action = na.exclude) #mod1_all has collinear variables

# ordiR2step
step_r2_all = ordiR2step(mod0_all, scope = formula(mod1_all), direction = "both", perm.max = 200, na.action = na.exclude) #mod1_all has collinear variables
print(step_r2_all)
# + host_species AIC = 2894.7, F = 1.783, p = 0.002

############ PREVIOUS CODE #############


# _ mammals ----
# by species 
rda_sp_mammals = capscale(formula = df_otus_sub_mammals ~ host_species, data = df_metadata_sub_mammals,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_sp_mammals = anova(rda_sp_mammals)
print(aov_rda_sp_mammals) #p = 0.001
RsquareAdj(rda_sp_mammals) #adj R^2 = 0.07303849

# by diet
rda_diet_mammals = capscale(formula = df_otus_sub_mammals ~ basic_diet, data = df_metadata_sub_mammals,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_diet_mammals = anova(rda_diet_mammals)
print(aov_rda_diet_mammals) #p = 0.001 
RsquareAdj(rda_diet_mammals) #adj R^2 = 0.008275982

# by sociality
rda_social_mammals = capscale(formula = df_otus_sub_mammals ~ basic_sociality, data = df_metadata_sub_mammals,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_social_mammals = anova(rda_social_mammals)
print(aov_rda_social_mammals) #p = 0.001 
RsquareAdj(rda_social_mammals) #adj R^2 = 0.02098

# model selection
mod0_mammals = capscale(df_otus_sub_mammals ~ 1, data = df_metadata_sub_mammals, distance = "robust.aitchison", na.action = na.exclude)
mod1_mammals = capscale(formula = df_otus_sub_mammals ~ host_species + basic_diet + basic_sociality, data = df_metadata_sub_mammals,  distance = "robust.aitchison", na.action = na.exclude) #some collinearity

# ordistep
step_r2_mammals = ordiR2step(mod0_mammals, scope = formula(mod1_mammals), perm.max = 200, na.action = na.exclude)
print(step_r2_mammals)

# _ birds ----
# by species
rda_sp_birds = capscale(formula = df_otus_sub_birds ~ host_species, data = df_metadata_sub_birds,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_sp_birds = anova(rda_sp_birds)
print(aov_rda_sp_birds) #p = 0.001
RsquareAdj(rda_sp_birds) #ajd R^2 = 0.1229736

# by diet
rda_diet_birds = capscale(formula = df_otus_sub_birds ~ basic_diet, data = df_metadata_sub_birds,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_diet_birds = anova(rda_diet_birds)
print(aov_rda_diet_birds) #p = 0.001
RsquareAdj(rda_diet_birds) #adj R^2 = 0.0224344

# by sociality
rda_social_birds = capscale(formula = df_otus_sub_birds ~ basic_sociality, data = df_metadata_sub_birds,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_social_birds = anova(rda_social_birds)
print(aov_rda_social_birds) #p = 0.012
RsquareAdj(rda_social_birds) #adj R^2 = 0.02076328

# _ model selection
# run null and full models
mod0_birds = capscale(df_otus_sub_birds ~ 1, data = df_metadata_sub_birds, distance = "robust.aitchison", na.action = na.exclude)
mod1_birds = capscale(formula = df_otus_sub_birds ~ host_species + basic_diet + basic_sociality, data = df_metadata_sub_birds,  distance = "robust.aitchison", na.action = na.exclude)

# ordistep
step_r2_birds = ordiR2step(mod0_birds, scope = formula(mod1_birds), perm.max = 200, na.action = na.exclude)
print(step_r2_birds)

# _ plot all hosts ----
df_metadata_sub = df_metadata_sub |> rownames_to_column() |> rename("sample_id" = "rowname")

# species
rda_scores_sp_df = rda_sp_all |> 
  scores(display = "sites") |> # extract the site scores
  as.data.frame() |> 
  rownames_to_column() |> 
  mutate(sample_id = rowname) |> #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_sp_all, groups = df_metadata_sub$host_species, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 4.35%, CAP2 1.39%

rda_scores_sp_df |> 
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
rda_scores_diet_df = rda_diet_all |> 
  scores(display = "sites") |> # extract the site scores
  as.data.frame() |> 
  rownames_to_column() |> 
  mutate(sample_id = rowname) |> #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_diet_all, groups = df_metadata_sub$basic_diet, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 1.13%, CAP2 0.67%

rda_scores_diet_df |> 
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
rda_scores_social_df = rda_social_all |> 
  scores(display = "sites") |> # extract the site scores
  as.data.frame() |> 
  rownames_to_column() |> 
  mutate(sample_id = rowname) |> #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_social_all, groups = df_metadata_sub$basic_sociality, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 1.66%, CAP2 0.52%

rda_scores_social_df |> 
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
rda_scores_sp_mammals_df = rda_sp_mammals |> 
  scores(display = "sites") |> # extract the site scores
  as.data.frame() |> 
  rownames_to_column() |> 
  mutate(sample_id = rowname) |> #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_sp_mammals, groups = df_metadata_sub_mammals$host_species, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 4.52%, CAP2 1.57%

rda_scores_sp_mammals_df |> 
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
rda_scores_diet_mammals_df = rda_diet_mammals |> 
  scores(display = "sites") |> # extract the site scores
  as.data.frame() |> 
  rownames_to_column() |> 
  mutate(sample_id = rowname) |> #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_diet_mammals, groups = df_metadata_sub_mammals$basic_diet, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 1.28%, CAP2 0.19%

rda_scores_diet_mammals_df |> 
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
rda_scores_social_mammals_df = rda_social_mammals |> 
  scores(display = "sites") |> # extract the site scores
  as.data.frame() |> 
  rownames_to_column() |> 
  mutate(sample_id = rowname) |> #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_social_mammals, groups = df_metadata_sub_mammals$basic_sociality, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 2.03%, CAP2 0.70%

rda_scores_social_mammals_df |> 
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
df_metadata_sub_birds = df_metadata_sub_birds |> rownames_to_column() |> rename("sample_id" = "rowname")

# species
rda_scores_sp_birds_df = rda_sp_birds |> 
  scores(display = "sites") |> # extract the site scores
  as.data.frame() |> 
  rownames_to_column() |> 
  mutate(sample_id = rowname) |> #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_sp_birds, groups = df_metadata_sub_birds$host_species, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 10.65%, CAP2 4.17%

rda_scores_sp_birds_df |> 
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
rda_scores_diet_birds_df = rda_diet_birds |> 
  scores(display = "sites") |> # extract the site scores
  as.data.frame() |> 
  rownames_to_column() |> 
  mutate(sample_id = rowname) |> #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_diet_birds, groups = df_metadata_sub_birds$basic_diet, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 3.82% (no other axes)

rda_scores_diet_birds_df |> 
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
rda_scores_social_birds_df = rda_social_birds |> 
  scores(display = "sites") |> # extract the site scores
  as.data.frame() |> 
  rownames_to_column() |> 
  mutate(sample_id = rowname) |> #add a sample_id column
  left_join(df_metadata_sub, by = "sample_id")

gg_ordiplot(rda_social_birds, groups = df_metadata_sub_birds$basic_sociality, hull = FALSE, label = FALSE, spiders = TRUE, ellipse = FALSE, pt.size = 2, plot = TRUE) #CAP1 4.07%, CAP2 1.16%

rda_scores_social_birds_df |> 
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


################ DIFFERENTIAL ABUNDANCE ################

##### setup #####

# all hosts
df_metadata_sub_da = df_metadata_sub |> 
  filter(basic_sociality != "intermediate") |> # exclude intermediate sociality
  column_to_rownames("sample_id") # bring rownames back

# subsample otu table
common_samples_sub_da = intersect(rownames(mx_otus), rownames(df_metadata_sub_da))

df_otus_sub_da = df_otus_sub |>
  rownames_to_column() |>
  filter(rowname %in% common_samples_sub_da) |> 
  column_to_rownames("rowname") |>
  select(where(~ sum(.) != 0))

# mammals
df_metadata_sub_mammals_da = df_metadata_sub |> 
  filter(basic_sociality != "intermediate") |>
  filter(host_class == "c__Mammalia") |> 
  column_to_rownames("sample_id")

common_samples_sub_mammals_da = intersect(rownames(mx_otus), rownames(df_metadata_sub_mammals_da))

df_otus_sub_mammals_da = df_otus_sub |>
  rownames_to_column() |>
  filter(rowname %in% common_samples_sub_mammals_da) |> 
  column_to_rownames("rowname") |>
  select(where(~ sum(.) != 0))

# birds
df_metadata_sub_birds_da = df_metadata_sub |> 
  filter(basic_sociality != "intermediate") |>
  filter(host_class == "c__Aves") |> 
  column_to_rownames("sample_id")

common_samples_sub_birds_da = intersect(rownames(mx_otus), rownames(df_metadata_sub_birds_da))

df_otus_sub_birds_da = df_otus_sub |>
  rownames_to_column() |>
  filter(rowname %in% common_samples_sub_birds_da) |> 
  column_to_rownames("rowname") |>
  select(where(~ sum(.) != 0))

# run this on the server only
# all hosts
v_social_da = df_metadata_sub_da |> pull(basic_sociality)
social_clr = aldex.clr(t(df_otus_sub_da), v_social_da, mc.samples = 200) #basic_sociality is a character vector, so aldex sorts it alphabetically -> social is the reference group
social_ttest = aldex.ttest(social_clr) #do not run on laptop
aldex_social_effect = aldex.effect(social_clr, CI = TRUE)
social_aldex_all = data.frame(social_ttest, aldex_social_effect)
social_aldex_all = social_aldex_all |> rownames_to_column()

# mammals
v_social_da_mammals = df_metadata_sub_mammals_da |> pull(basic_sociality)
social_clr_mammals = aldex.clr(t(df_otus_sub_mammals_da), v_social_da_mammals, mc.samples = 200)
social_ttest_mammals = aldex.ttest(social_clr_mammals)
aldex_social_effect_mammals = aldex.effect(social_clr_mammals, CI = TRUE)
social_aldex_mammals = data.frame(social_ttest_mammals, aldex_social_effect_mammals)
social_aldex_mammals = social_aldex_mammals |> rownames_to_column()

# birds
v_social_da_birds = df_metadata_sub_birds_da |> pull(basic_sociality)
social_clr_birds = aldex.clr(t(df_otus_sub_birds_da), v_social_da_birds, mc.samples = 200)
social_ttest_birds = aldex.ttest(social_clr_birds)
aldex_social_effect_birds = aldex.effect(social_clr_birds, CI = TRUE)
social_aldex_birds = data.frame(social_ttest_birds, aldex_social_effect_birds)
social_aldex_birds = social_aldex_birds |> rownames_to_column()

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
taxonomy_table = p_tax |> as.data.frame() |> rownames_to_column()

# generate df of differentially abundant OTUs
social_diff_asvs_mammals = social_aldex_mammals |> 
  filter(we.eBH < 0.001) |> 
  rownames_to_column("...1") |> 
  left_join(taxonomy_table, by = "rowname") |> # add ASV taxonomical info
  replace_na(list(genus = "")) |>
  replace_na(list(species = "sp.")) |>
  mutate(across(c("genus"), substr, 6, nchar(genus))) |> 
  mutate(across(c("species"), substr, 6, nchar(species))) |> 
  mutate(ci_product = effect.low * effect.high ) |> #column that returns TRUE if CI does not contain 0# |> 
  mutate(ci_no_zero = ci_product > 0) #returns TRUE if CI does not include zero

social_diff_asvs_mammals |> 
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
social_diff_asvs_birds = social_aldex_birds |>
  filter(we.eBH < 0.05) |> #we.eBH = expected Benjamini-Hochberg corrected p-value of Welch’s t test (higher p still yields 0 differentially abundant OTUs)
  rownames_to_column("...1") |> 
  left_join(taxonomy_table, by = "rowname") |> # add ASV taxonomical info
  replace_na(list(genus = "")) |>
  replace_na(list(species = "sp.")) |>
  mutate(otu_scientific_name = paste(genus, species, sep = " ")) |> #concatenate to full ASV name
  mutate(otu_scientific_name = str_trim(otu_scientific_name, side = "left")) |>
  mutate(otu_scientific_name = as.factor(otu_scientific_name)) |> 
  mutate(ci_product = effect.low * effect.high ) |> #column that returns TRUE if CI does not contain 0# |> 
  mutate(ci_no_zero = ci_product > 0) #returns TRUE if CI does not include zero
#there are no differentially abundant OTUs in birds