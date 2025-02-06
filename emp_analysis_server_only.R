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
library(jtools)
library(vegan)
library(lme4)
library(ALDEx2)
theme_set(theme_classic())
filter = dplyr::filter
select = dplyr::select

## file input ----
load(file = "emp_analysis_data_for_server.RData")

# data files
#df_metadata_sub, df_otus_sub, ps_sub, df_metadata_sub_mammals, df_metadata_sub_birds, df_metadata_sub_birds, df_otus_sub_birds, ps_mammals_sub, ps_birds_sub, dissimilarity_table_mammals, dissimilarity_table_birds

## relative abundance ---- 
# mammals
total_mammals = median(sample_sums(ps_mammals_sub))
standf = function(x, t = total_mammals) round(t * (x / sum(x)))
ps_norm_mammals = transform_sample_counts(ps_mammals_sub, standf)

# convert ps to a tidy data frame
relab_mammals_df = ps_norm_mammals |>
  psmelt() |>
  group_by(Sample, phylum) |>
  summarize(Abundance = sum(Abundance), .groups = "drop") |>
  mutate(RelativeAbundance = Abundance / sum(Abundance))

# get top 10 phyla
top_phyla_mammals = relab_mammals_df |>
  group_by(phylum) |>
  summarize(TotalAbundance = sum(Abundance), .groups = "drop") |>
  arrange(desc(TotalAbundance)) |>
  slice_head(n = 10) |>
  pull(phylum)

# print top 10 phyla
print(top_phyla_mammals)

# reclassify other phyla as "Other"
relab_mammals_df = relab_mammals_df |>
  mutate(phylum = ifelse(phylum %in% top_phyla_mammals, phylum, "Other"))

# recalculate relative abundances with new categories
relab_mammals_df = relab_mammals_df |>
  group_by(Sample, phylum) |>
  summarize(Abundance = sum(Abundance), .groups = "drop") |>
  mutate(RelativeAbundance = Abundance / sum(Abundance))

# plot
plot_relab_mammals = relab_mammals_plot = ggplot(relab_mammals_df, aes(x = Sample, y = RelativeAbundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab("sample") +
  ylab("relative abundance") +
  theme(axis.text.x = element_blank(),
        text = element_text(size = 14)) +
  scale_fill_viridis_d(labels = c("Proteobacteria", "Firmicutes", "Bacteroidetes", "Actinobacteria", "Actinobacteria", "Tenericutes", "Verrucomicrobia", "Lentisphaerae", "Spirochaetae", "Euryarchaeota", "Cyanobacteria"))
plot_relab_mammals

save(relab_mammals_df, top_phyla_mammals, plot_relab_mammals, relab_mammals_df, file = "plot_relab_mammals.RData")

## dissimilarity stats ----
glm_diss_mammals = lmer(diss ~ basic_sociality + (1|host_species_1), data = dissimilarity_table_mammals)
print("glm_diss_mammals summary----")
summary(glm_diss_mammals)
print("glm_diss_mammals jtools----")
jtools::summ(glm_diss_mammals) #Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()

## diff ab ----
# all
df_metadata_sub_da = df_metadata_sub |> 
  filter(basic_sociality != "intermediate") #|> # exclude intermediate sociality
  #column_to_rownames("sample_id") # bring rownames back

# subsample otu table
common_samples_sub_da = intersect(rownames(df_otus_sub), rownames(df_metadata_sub_da))

df_otus_sub_da = df_otus_sub |>
  rownames_to_column() |>
  filter(rowname %in% common_samples_sub_da) |> 
  column_to_rownames("rowname") |>
  select(where(~ sum(.) != 0))

# mammals
df_metadata_sub_mammals_da = df_metadata_sub |> 
  filter(basic_sociality != "intermediate") |>
  filter(host_class == "c__Mammalia") #|> 
  #column_to_rownames("sample_id")

common_samples_sub_mammals_da = intersect(rownames(df_otus_sub), rownames(df_metadata_sub_mammals_da))

df_otus_sub_mammals_da = df_otus_sub |>
  rownames_to_column() |>
  filter(rowname %in% common_samples_sub_mammals_da) |> 
  column_to_rownames("rowname") |>
  select(where(~ sum(.) != 0))

# birds
df_metadata_sub_birds_da = df_metadata_sub |> 
  filter(basic_sociality != "intermediate") |>
  filter(host_class == "c__Aves") #|> 
  #column_to_rownames("sample_id")

common_samples_sub_birds_da = intersect(rownames(df_otus_sub), rownames(df_metadata_sub_birds_da))

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
social_aldex_all = social_aldex_all |> rownames_to_column()

# _ mammals
v_social_da_mammals = df_metadata_sub_mammals_da |> pull(basic_sociality)
social_clr_mammals = aldex.clr(t(df_otus_sub_mammals_da), v_social_da_mammals, mc.samples = 200)
social_ttest_mammals = aldex.ttest(social_clr_mammals)
aldex_social_effect_mammals = aldex.effect(social_clr_mammals, CI = TRUE)
social_aldex_mammals = data.frame(social_ttest_mammals, aldex_social_effect_mammals)
social_aldex_mammals = social_aldex_mammals |> rownames_to_column()

# _ birds
v_social_da_birds = df_metadata_sub_birds_da |> pull(basic_sociality)
social_clr_birds = aldex.clr(t(df_otus_sub_birds_da), v_social_da_birds, mc.samples = 200)
social_ttest_birds = aldex.ttest(social_clr_birds)
aldex_social_effect_birds = aldex.effect(social_clr_birds, CI = TRUE)
social_aldex_birds = data.frame(social_ttest_birds, aldex_social_effect_birds)
social_aldex_birds = social_aldex_birds |> rownames_to_column()

write_csv(social_aldex_all, "social_aldex_all.csv")
write_csv(social_aldex_mammals, "social_aldex_mammals.csv")
write_csv(social_aldex_birds, "social_aldex_birds.csv")

save.image("emp_analysis_server_only_results.RData")