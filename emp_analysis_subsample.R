#### 
#
# emp analysis on subsampled data
#
####

## packages ----
library(tidyverse)
library(conflicted)
library(viridis)
library(phyloseq)
library(fantaxtic)
library(ggpubr)
library(patchwork)
library(vegan)
library(lme4)
library(ALDEx2)
theme_set(theme_classic())
filter = dplyr::filter
select = dplyr::select

## file input ----
load(file = "emp_analysis_data_2025_02_03.RData") #this dataset has been modified to omit 1 reptile species, black bear, colobine primates

## subsampling ----
# _ dfs ----
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
df_metadata_sub <- df_metadata_sub |> 
  slice(match(rownames(df_otus_sub), rownames(.)))
all(rownames(df_otus_sub) == rownames(df_metadata_sub)) #TRUE

# _ phyloseq ----
# set up metadata & otu table
p_metadata_sub = sample_data(df_metadata_sub)
p_otu_sub = otu_table(df_otus_sub, taxa_are_rows = FALSE)

# check that sample names match
sample_names(p_metadata_sub)
sample_names(p_otu_sub)

# assemble phyloseq object
physeq_sub = phyloseq(p_otu_sub, p_metadata_sub, p_tax)
ps_sub = merge_phyloseq(physeq_sub, tree_emp)

## subset classes ----
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

## descriptive stats ----
# reads
# total number of reads
df_otus_sub |> sum() #12,950,659 (this number will vary slightly each time this script is run since subsampling is random)

# number of reads per sample
mean(rowSums(df_otus_sub)) #34,535.09
sd(rowSums(df_otus_sub)) #23,192.04

# _ all ----
df_metadata_sub |> nrow() #375 samples
df_metadata_sub |> group_by(host_species) |> 
  summarise(n = n()) |> print(n = 100) #across 45 species
df_otus_sub |> as.data.frame() |> ncol() #34,482

# _ mammals ----
df_metadata_sub_mammals |> nrow() #312 samples
df_metadata_sub_mammals |> group_by(host_species) |> 
  summarise(n = n()) |> print(n = 35) #across 35 species
df_otus_sub_mammals |> ncol() #31,451 OTUs

# _ birds ----
df_metadata_sub_birds |> nrow() #63 samples
df_metadata_sub_birds |> group_by(host_species) |> 
  summarise(n = n()) |> print(n = 100) #across 10 species
df_otus_sub_birds |> ncol() #13,915

## relative abundances ---- 
# _ SERVER plots ----
# normalize number of reads using median sequencing depth
total = median(sample_sums(ps_sub))
standf = function(x, t = total) round(t * (x / sum(x)))
ps_norm = transform_sample_counts(ps_sub, standf)

# plot relative abundances per class (run mammals on server)
# mammals - does not run on laptop
total_mammals = median(sample_sums(ps_mammals_sub))
standf = function(x, t = total_mammals) round(t * (x / sum(x)))
ps_norm_mammals = transform_sample_counts(ps_mammals_sub, standf)

relab_mammals = plot_bar(ps_norm_mammals, fill = "phylum") +
  geom_bar(aes(colour = phylum, fill = phylum), stat = "identity", position = "stack") +
  theme(legend.position = "none") #+ facet_grid(~ host_order)

# birds
total_birds = median(sample_sums(ps_birds_sub))
standf = function(x, t = total_birds) round(t * (x / sum(x)))
ps_norm_birds = transform_sample_counts(ps_birds_sub, standf)

relab_birds = plot_bar(ps_norm_birds, fill = "phylum") +
  geom_bar(aes(colour = phylum, fill = phylum), stat = "identity", position = "stack") +
  theme(legend.position = "none") #+ facet_grid(~ host_order)

# _ dominant taxa ----
# mammals
fantaxtic::top_taxa(ps_mammals_sub, n = 3, tax_level = "phylum")
# Firmicutes 0.367, Proteobacteria 0.334, Bacteroidetes 0.180
fantaxtic::top_taxa(ps_mammals_sub, n = 5, tax_level = "class")

# birds
fantaxtic::top_taxa(ps_birds_sub, n = 3, tax_level = "phylum")
#Firmicutes 0.396, Proteobacteria 0.338, Bacteroidetes 0.144)
fantaxtic::top_taxa(ps_birds_sub, n = 5, tax_level = "class")
#Bacilli 0.211, Clostridia 0.175, Gammaproteobacteria 0.173, Betaproteobacteria 0.0645, Alphaproteobacteria 0.0730

## alpha div ----
# _ class ----
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

# _ sociality ----
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

kruskal.test(adiv_faith_pd ~ basic_sociality, data = df_metadata_sub) #KW chi-squared = 15.286, df = 2, p-value = 0.0004795
kruskal.test(adiv_observed_otus ~ basic_sociality, data = df_metadata_sub) #KW chi-squared = 12.594, df = 2, p-value = 0.001842

# sociality is not independent of host phylogeny
glm_soc_order = lmer(adiv_faith ~ basic_sociality + (1|host_order), data = df_metadata_sub)
summary(glm_soc_order)
jtools::summ(glm_soc_order) #not significant for faith, obs_otus, chao1, shannon

## dissimilarity ----
# _ mammals ----
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

# _ birds ----
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

# _ plots ----
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

# _ SERVER stats ----
glm_diss_mammals = lmer(diss ~ basic_sociality + (1|host_species_1), data = dissimilarity_table_mammals)
summary(glm_diss_mammals)
jtools::summ(glm_diss_mammals) #Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()

glm_diss_birds = lmer(diss ~ basic_sociality + (1|host_species_1), data = dissimilarity_table_birds)
summary(glm_diss_birds)
jtools::summ(glm_diss_birds) #social p = 0.60, solitary 0.46


# db-RDA ----
# _ all hosts ----
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
print(aov_rda_diet_all) #p = 0.001
RsquareAdj(rda_social_all)#ajd R^2 = 0.01653155

# _ model selection ----
# run null and full models
mod0_all = capscale(df_otus_sub ~ 1, data = df_metadata_sub, na.action = na.exclude)
mod1_all = capscale(formula = df_otus_sub ~ host_species + basic_diet + basic_sociality, data = df_metadata_sub,  distance = "robust.aitchison", na.action = na.exclude) #Some constraints or conditions were aliased because they were redundant. This can happen if terms are linearly dependent (collinear): ‘host_speciess__Turdus_olivater’, ‘host_speciess__Vulpes_vulpes’, ‘host_speciess__Zonotrichia_capensis’

# ordistep
step_r2_all = ordiR2step(mod0_all, scope = formula(mod1_all), perm.max = 200, na.action = na.exclude) #mod1_all has colinear variables
print(step_r2_all)

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

# model selection ----
mod0_mammals = capscale(df_otus_sub_mammals ~ 1, data = df_metadata_sub_mammals, na.action = na.exclude)
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

# _ model selection ----
# run null and full models
mod0_birds = capscale(df_otus_sub_birds ~ 1, data = df_metadata_sub_birds, na.action = na.exclude)
mod1_birds = capscale(formula = df_otus_sub_birds ~ host_species + basic_diet + basic_sociality, data = df_metadata_sub_birds,  distance = "robust.aitchison", na.action = na.exclude)

# ordistep
step_r2_birds = ordiR2step(mod0_birds, scope = formula(mod1_birds), perm.max = 200, na.action = na.exclude)
print(step_r2_birds)

# plot rda ----
# _ all hosts ----
df_metadata_sub = df_metadata_sub |> rownames_to_column() |> rename("sample_id" = "rowname")

# species
rda_scores_sp = scores(rda_sp_all, display = "sites") # extract the site scores
rda_scores_sp_df = as.data.frame(rda_scores_sp)
rda_scores_sp_df$sample_id = rownames(rda_scores_sp_df) # add a study_id column
rda_fort_sp_all = rda_scores_sp_df |> 
  left_join(df_metadata_sub, by = "sample_id") # join to metadata

# diet
rda_scores_diet = scores(rda_diet_all, display = "sites")
rda_scores_diet_df = as.data.frame(rda_scores_diet)
rda_scores_diet_df$sample_id = rownames(rda_scores_diet_df)
rda_fort_diet_all = rda_scores_diet_df |> 
  left_join(df_metadata_sub, by = "sample_id")

# sociality
rda_scores_social = scores(rda_social_all, display = "sites")
rda_scores_social_df = as.data.frame(rda_scores_social)
rda_scores_social_df$sample_id = rownames(rda_scores_social_df)
rda_fort_social_all = rda_scores_social_df |> 
  left_join(df_metadata_sub, by = "sample_id")

# TODO ----
# plots
rda_fort_sp_all |> 
  ggplot(aes(x = CAP1, y = CAP2, colour = host_species)) +
  geom_point() +
  theme(text = element_text(size = 14))

rda_fort_diet_all |> 
  ggplot(aes(x = CAP1, y = CAP2, colour = basic_diet)) +
  geom_point() +
  theme(text = element_text(size = 14))

rda_fort_social_all |> 
  ggplot(aes(x = CAP1, y = CAP2, colour = host_class)) +
  geom_point() +
  theme(text = element_text(size = 14))

p_sp_all
p_diet_all
p_social_all

# _ mammals ----
df_metadata_sub_mammals = df_metadata_sub_mammals |> rownames_to_column() |> rename("sample_id" = "rowname")

# species
rda_scores_sp_mammals = scores(rda_sp_mammals, display = "sites")
rda_scores_sp_mammals_df = as.data.frame(rda_scores_sp_mammals)
rda_scores_sp_mammals_df$sample_id = rownames(rda_scores_sp_mammals_df)
rda_fort_sp_mammals_all = rda_scores_sp_mammals_df |> 
  left_join(df_metadata_sub_mammals, by = "sample_id")

# diet
rda_scores_diet_mammals = scores(rda_diet_mammals, display = "sites")
rda_scores_diet_mammals_df = as.data.frame(rda_scores_diet_mammals)
rda_scores_diet_mammals_df$sample_id = rownames(rda_scores_diet_mammals_df)
rda_fort_diet_mammals_all = rda_scores_diet_mammals_df |> 
  left_join(df_metadata_sub_mammals, by = "sample_id")

# sociality
rda_scores_social_mammals = scores(rda_social_mammals, display = "sites")
rda_scores_social_mammals_df = as.data.frame(rda_scores_social_mammals)
rda_scores_social_mammals_df$sample_id = rownames(rda_scores_social_mammals_df)
rda_fort_social_mammals = rda_scores_social_mammals_df |> 
  left_join(df_metadata_sub_mammals, by = "sample_id")

# plots
p_sp_mammals = ggplot(rda_fort_social_mammals, aes(x = CAP1, y = CAP2, colour = host_species)) +
  geom_point()
p_diet_mammals = ggplot(rda_fort_social_mammals, aes(x = CAP1, y = CAP2, colour = basic_diet)) +
  geom_point()
p_social_mammals = ggplot(rda_fort_social_mammals, aes(x = CAP1, y = CAP2, colour = host_species)) +
  geom_point() + theme(legend.position="none")

# _ birds ----
df_metadata_sub_birds = df_metadata_sub_birds |> rownames_to_column() |> rename("sample_id" = "rowname")

# species
rda_scores_sp_birds = scores(rda_sp_birds, display = "sites") # extract the site scores
rda_scores_sp_df_birds = as.data.frame(rda_scores_sp_birds)
rda_scores_sp_df_birds$sample_id = rownames(rda_scores_sp_df_birds) # add a study_id column
rda_fort_sp_birds = rda_scores_sp_df_birds |> 
  left_join(df_metadata_sub_birds, by = "sample_id") # join to metadata

# diet
rda_scores_diet_birds = scores(rda_diet_birds, display = "sites")
rda_scores_diet_df_birds = as.data.frame(rda_scores_diet_birds)
rda_scores_diet_df_birds$sample_id = rownames(rda_scores_diet_df_birds)
rda_fort_diet_birds = rda_scores_diet_df_birds |> 
  left_join(df_metadata_sub, by = "sample_id")

# sociality
rda_scores_social_birds = scores(rda_social_birds, display = "sites")
rda_scores_social_df_birds = as.data.frame(rda_scores_social_birds)
rda_scores_social_df_birds$sample_id = rownames(rda_scores_social_df_birds)
rda_fort_social_birds = rda_scores_social_df_birds |> 
  left_join(df_metadata_sub_birds, by = "sample_id")

# plots
p_sp_birds = ggplot(rda_fort_sp_birds, aes(x = CAP1, y = CAP2, colour = host_species)) +
  geom_point() #+ theme(legend.position="none")
p_diet_birds = ggplot(rda_fort_diet_birds, aes(x = CAP1, y = CAP2, colour = basic_diet)) +
  geom_point()
p_social_birds = ggplot(rda_fort_social_birds, aes(x = CAP1, y = CAP2, colour = basic_sociality)) +
  geom_point()


# diff ab ----
# _ setup ----
# exclude intermediate sociality
df_metadata_sub_da = df_metadata_sub |> 
  filter(basic_sociality != "intermediate")
df_metadata_sub_mammals_da = df_metadata_sub_mammals |> 
  filter(basic_sociality != "intermediate")
df_metadata_sub_birds_da = df_metadata_sub_birds |> 
  filter(basic_sociality != "intermediate")

# subset otu table
common_samples_sub_da = intersect(rownames(df_otus_sub), rownames(df_metadata_sub_da))
df_otus_sub_da = df_otus_sub[common_samples_sub_da, , drop = FALSE]
df_otus_sub_da = df_otus_sub_da |> as.data.frame() |> select(where(~ sum(.) != 0)) |> #remove columns (OTUs that are 0 everywhere)
  as.matrix()

common_samples_sub_mammals_da = intersect(rownames(df_otus_sub_mammals), rownames(df_metadata_sub_mammals_da))
df_otus_sub_mammals_da = df_otus_sub_mammals[common_samples_sub_mammals_da, , drop = FALSE]
df_otus_sub_mammals_da = df_otus_sub_mammals_da |> as.data.frame() |> select(where(~ sum(.) != 0)) |> #remove columns (OTUs that are 0 everywhere)
  as.matrix()

common_samples_sub_birds_da = intersect(rownames(df_otus_sub_birds), rownames(df_metadata_sub_birds_da))
df_otus_sub_birds_da = df_otus_sub_birds[common_samples_sub_birds_da, , drop = FALSE]
df_otus_sub_birds_da = df_otus_sub_birds_da |> as.data.frame() |> select(where(~ sum(.) != 0)) |> #remove columns (OTUs that are 0 everywhere)
  as.matrix()

# _ all hosts ----
v_social_da = df_metadata_sub |> pull(basic_sociality)
social_clr = aldex.clr(t(df_otus_sub_da), v_social_da, mc.samples = 200) #do not run on laptop
social_ttest = aldex.ttest(social_clr) #do not run on laptop
aldex_social_effect = aldex.effect(social_clr, CI = TRUE)
social_aldex_all = data.frame(social_ttest, aldex_social_effect)

# _ mammals ----
v_social_da_mammals = df_metadata_sub_mammals_da |> pull(basic_sociality)
social_clr_mammals = aldex.clr(t(df_otus_sub_mammals_da), v_social_da_mammals, mc.samples = 200)
social_ttest_mammals = aldex.ttest(social_clr_mammals)
aldex_social_effect_mammals = aldex.effect(social_clr_mammals, CI = TRUE)
social_aldex_mammals = data.frame(social_ttest_mammals, aldex_social_effect_mammals)

# _ birds ----
v_social_da_birds = df_metadata_sub_birds_da |> pull(basic_sociality)
social_clr_birds = aldex.clr(t(df_otus_sub_birds_da), v_social_da_birds, mc.samples = 200)
social_ttest_birds = aldex.ttest(social_clr_birds)
aldex_social_effect_birds = aldex.effect(social_clr_birds, CI = TRUE)
social_aldex_birds = data.frame(social_ttest_birds, aldex_social_effect_birds)

write_csv(social_aldex_all, "social_aldex_all.csv")
write_csv(social_aldex_mammals, "social_aldex_mammals.csv")
write_csv(social_aldex_birds, "social_aldex_birds.csv")

save(social_aldex_all, social_aldex_mammals, social_aldex_birds, file = "diff_ab_data.RData")


# _ plots ----
par(mfrow = c(1,2))
aldex.plot(social_aldex_all, type = "MA", test = "welch", main = "MA plot")
aldex.plot(social_aldex_all, type = "MW", test = "welch", main = "effect plot")

# identify OTUs by name
taxonomy_table = p_tax |> as.data.frame() |> rownames_to_column()

# generate df of differentially abundant ASVs
# mammals only, since there are none in birds (but keep the code for birds below)
social_diff_asvs_mammals = social_aldex_mammals |> as.data.frame() |> rownames_to_column() |> 
  filter(we.eBH < 0.05) |> 
  rownames_to_column("...1") |> 
  left_join(taxonomy_table, by = "rowname") |> # add ASV taxonomical info
  replace_na(list(genus = "")) |>
  replace_na(list(species = "sp.")) |>
  mutate(otu_scientific_name = paste(genus, species, sep = " ")) |> #concatenate to full ASV name
  mutate(otu_scientific_name = str_trim(otu_scientific_name, side = "left")) |>
  mutate(otu_scientific_name = as.factor(otu_scientific_name)) |> 
  mutate(ci_product = effect.low * effect.high ) |> #column that returns TRUE if CI does not contain 0# |> 
  mutate(ci_no_zero = ci_product > 0) #returns TRUE if CI does not include zero

# birds
social_diff_asvs_birds = social_aldex_birds |> as.data.frame() |> rownames_to_column() |> 
  filter(we.eBH < 0.05) |> 
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

social_diff_asvs_mammals |> 
  ggplot(aes(x = effect, y = reorder(otu_scientific_name, (effect*-1)), colour = phylum)) +
  geom_point() +
  xlab("Effect size") +
  ylab("ASV") +
  theme(axis.text.y = element_text(face = "italic")) +
  ggtitle("Differentially abundant ASVs in mammals\nSocial and solitary\n(BH-corrected)") + 
  geom_vline(xintercept = 0, linetype = "dotted")

