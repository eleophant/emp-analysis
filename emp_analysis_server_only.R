#### 
#
# emp analysis on subsampled data - server-run code
#
####

## packages ----
library(tidyverse)
library(conflicted)
library(viridis)
library(phyloseq)
library(ggordiplots)
library(vegan)
library(lme4)
library(ALDEx2)
theme_set(theme_classic())
filter = dplyr::filter
select = dplyr::select

## file input ----
load(file = "emp_data_for_server.RData") #this dataset has been modified to omit 1 reptile species, black bear, colobine primates

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
df_metadata_sub = df_metadata_sub |>
  slice(match(rownames(df_otus_sub), rownames(df_metadata_sub)))
all(rownames(df_otus_sub) == rownames(df_metadata_sub)) #TRUE

# _ phyloseq ----
# set up metadata & otu table
p_metadata_sub = sample_data(df_metadata_sub)
p_otu_sub = otu_table(df_otus_sub, taxa_are_rows = FALSE)

# check that sample names match
all(sample_names(p_metadata_sub) == sample_names(p_otu_sub)) #TRUE

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

## relative abundance ---- 
# _ plots ----
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

library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats) # For better categorical handling

# Convert phyloseq object to a tidy data frame
df_relab_birds <- ps_birds_sub %>%
  phyloseq::psmelt() %>%  # Converts phyloseq object to a tidy format
  group_by(Sample, phylum) %>%  
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%  # Summarize by phylum
  mutate(Relative_Abundance = Abundance / sum(Abundance))  # Convert to relative abundance

ggplot(df_relab_birds, aes(x = Sample, y = Relative_Abundance, fill = fct_reorder(phylum, -Relative_Abundance))) +
  geom_bar(stat = "identity", position = "stack") +  
  scale_fill_viridis_d() +
  theme_minimal() +  
  theme(axis.text.x = element_blank(),
        legend.position = "none") +  
  labs(x = "Sample", y = "Relative Abundance", fill = "Phylum")  

## pair dissimilarity ----
glm_diss_mammals = lmer(diss ~ basic_sociality + (1|host_species_1), data = dissimilarity_table_mammals)
summary(glm_diss_mammals)
jtools::summ(glm_diss_mammals) #Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()

glm_diss_birds = lmer(diss ~ basic_sociality + (1|host_species_1), data = dissimilarity_table_birds)
summary(glm_diss_birds)
jtools::summ(glm_diss_birds) #social p = 0.60, solitary 0.46

## diff ab ----
# _ setup ----
# all
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

# _ all hosts ----
v_social_da = df_metadata_sub_da |> pull(basic_sociality)
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
