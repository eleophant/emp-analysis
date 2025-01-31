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
#library(phyloseqCompanion)
library(vegan)
library(lme4)
library(ALDEx2)
theme_set(theme_classic())
filter <- dplyr::filter

# file input ----
load(file = "emp_analysis_rda_2024_06_18.RData")
df_otus = df_otus |> as.data.frame()

# subsampling ----
# _ dfs ----
df_obs_n = df_metadata |> 
  group_by(host_species) |> 
  summarise(n = n()) |> 
  print(n = 27)

df_obs_n |> summarise(mean(n), median(n), max(n), min(n))
#mean 41, median 3, max 902, min 1

# subsample
df_metadata_sub = df_metadata |> 
  rownames_to_column() |> 
  group_by(host_species) |> 
  slice_sample(n = 41)  |>
  ungroup() |> 
  column_to_rownames("rowname")

common_samples_sub = intersect(rownames(df_otus), rownames(df_metadata_sub))
df_otus_sub = df_otus[common_samples_sub, , drop = FALSE]
df_otus_sub = df_otus_sub |> as.data.frame() |> select(where(~ sum(.) != 0)) |> #remove columns (OTUs that are 0 everywhere)
  as.matrix()

# _ phyloseq ----
#set up metadata & otu table
p_metadata_sub = sample_data(df_metadata_sub)
p_otu_sub = otu_table(df_otus_sub, taxa_are_rows = FALSE)

# check that sample names match
sample_names(p_metadata_sub)
sample_names(p_otu_sub)

# assemble phyloseq object
physeq_sub = phyloseq(p_otu_sub, p_metadata_sub, p_tax) #works!
ps_sub = merge_phyloseq(physeq_sub, tree_emp)

# subset classes ----
# mammals
df_metadata_sub_mammals = df_metadata_sub |>
  filter(host_class == "c__Mammalia")
common_samples_sub_mammals = intersect(rownames(df_otus), rownames(df_metadata_sub_mammals))
df_otus_sub_mammals = df_otus_mammals[common_samples_sub_mammals, , drop = FALSE]
df_otus_sub_mammals = df_otus_sub_mammals |> as.data.frame() |> select(where(~ sum(.) != 0)) |> #remove empty columns (OTUs that are 0 everywhere)
  as.matrix()

# birds
df_metadata_sub_birds = df_metadata_sub |>
  filter(host_class == "c__Aves")
common_samples_sub_birds = intersect(rownames(df_otus), rownames(df_metadata_sub_birds))
df_otus_sub_birds = df_otus_sub_birds[common_samples_sub_birds, , drop = FALSE]
df_otus_sub_birds = df_otus_sub_birds |> as.data.frame() |> select(where(~ sum(.) != 0)) |> #remove empty columns (OTUs that are 0 everywhere)
  as.matrix()

# phyloseq
ps_mammals_sub = subset_samples(ps_sub, host_class == "c__Mammalia")
ps_birds_sub = subset_samples(ps_sub, host_class == "c__Aves")

# descriptive stats ----
# total number of reads
df_otus_sub |> sum() #12,950,659

# number of reads per sample
mean(rowSums(df_otus_sub)) #34,535.09
sd(rowSums(df_otus_sub)) #23,192.04

# number of ASVs (number of columns)
df_otus_sub |> as.data.frame() |> ncol() #34,482


## TODO relative abundances ---- 
# normalize number of reads using median sequencing depth
total = median(sample_sums(ps))
standf = function(x, t = total) round(t * (x / sum(x)))
a_ps_norm = transform_sample_counts(ps, standf)

# plot relative abundances *** does not run on laptop***
#plot_bar(a_ps_norm, fill = "phylum") +
#  geom_bar(aes(colour = phylum, fill = phylum), stat = "identity", position = "stack") #+ facet_grid(~ host_order.x)

# _ dominant taxa ----
fantaxtic::top_taxa(ps_birds, n = 3, tax_level = "phylum")
#Firmicutes (39.6%), Proteobacteria (33.8%), and Bacteroidetes (14.4%)
fantaxtic::top_taxa(ps, n = 5, tax_level = "class")
#Clostridia (24.8%), Gammaproteobacteria (12.4%), and Bacilli (10.6%)

# mammals
fantaxtic::top_taxa(ps_mammals, n = 3, tax_level = "phylum")
fantaxtic::top_taxa(ps_mammals, n = 5, tax_level = "class")

# birds
fantaxtic::top_taxa(ps_birds, n = 3, tax_level = "phylum")
fantaxtic::top_taxa(ps_birds, n = 5, tax_level = "class")

## TODO alpha div ----
# per host class
df_metadata |> 
  filter(host_class != "c__Reptilia") |> 
  ggplot(aes(x = host_class, y = adiv_faith_pd)) +
  geom_boxplot() +
  xlab("Host class") +
  ylab("Faith PD") +
  scale_x_discrete(labels=c("c__Aves" = "Aves", "c__Mammalia" = "Mammalia"))

# statistical tests
kruskal.test(df_metadata_mb$adiv_faith_pd, df_metadata_mb$host_class) #p = 1.108e-06
summary(glm(df_metadata_mb$adiv_chao1 ~ df_metadata_mb$host_class)) #p = 8.28e-07


# db-RDA ----
# _ all hosts ----
# by species
print("RDA species, all hosts ---------------")
rda_sp_all = capscale(formula = df_otus_sub ~ host_species, data = df_metadata_sub,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_sp_all = anova(rda_sp_all) 
print(aov_rda_sp_all) #p = 0.001
RsquareAdj(rda_sp_all) #R^2 = 0.02556273

# by diet
print("RDA diet, all hosts ---------------")
rda_diet_all = capscale(formula = df_otus_sub ~ basic_diet, data = df_metadata_sub,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_diet_all = anova(rda_diet_all)
print(aov_rda_diet_all) #p = 0.001
RsquareAdj(rda_diet_all) #R^2 = 0.005332206

# by sociality
print("RDA sociality, all hosts ---------------")
rda_social_all = capscale(formula = df_otus_sub ~ basic_sociality, data = df_metadata_sub,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_social_all = anova(rda_social_all)
print(aov_rda_diet_all) #p = 0.001
RsquareAdj(rda_social_all)  #R^2 = 0.0141168

# _ model selection ----
# run null and full models
mod0_all = capscale(df_otus_sub ~ 1, data = df_metadata_sub, na.action = na.exclude)
mod1_all = capscale(formula = df_otus_sub ~ species + diet + sociality, data = df_metadata_sub,  distance = "robust.aitchison", na.action = na.exclude)

# ordistep
step_r2_all = ordiR2step(mod0_all, scope = formula(mod1_all), perm.max = 200, na.action = na.exclude)
print(step_r2_all)

# _ mammals ----
# by species 
rda_sp_mammals = capscale(formula = df_otus_sub_mammals ~ host_species, data = df_metadata_sub_mammals,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_sp_mammals = anova(rda_sp_mammals)
print(aov_rda_sp_mammals) #p = 0.001
RsquareAdj(rda_sp_mammals) #R^2 = 0.02287834

# by diet
rda_diet_mammals = capscale(formula = df_otus_sub_mammals ~ basic_diet, data = df_metadata_sub_mammals,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_diet_mammals = anova(rda_diet_mammals)
print(aov_rda_diet_mammals) #p = 0.001 
RsquareAdj(rda_diet_mammals) #R^2 = 0.003841225

# by sociality
rda_social_mammals = capscale(formula = df_otus_sub_mammals ~ basic_sociality, data = df_metadata_sub_mammals,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_social_mammals = anova(rda_social_mammals)
print(aov_rda_social_mammals) #p = 0.001 
RsquareAdj(rda_social_mammals)  #R^2 = 0.0142302

# _ model selection ----
# run null and full models
mod0_mammals = capscale(df_otus_sub_mammals ~ 1, data = df_metadata_sub_mammals, na.action = na.exclude)
mod1_mammals = capscale(formula = df_otus_sub_mammals ~ species + diet + sociality, data = df_metadata_sub_mammals,  distance = "robust.aitchison", na.action = na.exclude)

# ordistep
step_r2_mammals = ordiR2step(mod0_mammals, scope = formula(mod1_mammals), perm.max = 200, na.action = na.exclude)
print(step_r2_mammals)


# _ birds ----
# by species
rda_sp_birds = capscale(formula = df_otus_sub_birds ~ host_species, data = df_metadata_sub_birds,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_sp_birds = anova(rda_sp_birds)
print(aov_rda_sp_birds) #p = 0.001
RsquareAdj(rda_sp_birds)  #R^2 = 0.1229736

# by diet
rda_diet_birds = capscale(formula = df_otus_sub_birds ~ basic_diet, data = df_metadata_sub_birds,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_diet_birds = anova(rda_diet_birds)
print(aov_rda_diet_birds) #p = 0.001
RsquareAdj(rda_diet_birds)  #R^2 = 0.0224344

# by sociality
rda_social_birds = capscale(formula = df_otus_sub_birds ~ basic_sociality, data = df_metadata_sub_birds,  distance = "robust.aitchison", na.action = na.exclude)
aov_rda_social_birds = anova(rda_social_birds)
print(aov_rda_social_birds) #p = 0.017
RsquareAdj(rda_social_birds)  #R^2 = 0.02076328


# _ model selection ----
# run null and full models
mod0_birds = capscale(df_otus_sub_birds ~ 1, data = df_metadata_sub_birds, na.action = na.exclude)
mod1_birds = capscale(formula = df_otus_sub_birds ~ species + diet + sociality, data = df_metadata_sub_birds,  distance = "robust.aitchison", na.action = na.exclude)

# ordistep
step_r2_birds = ordiR2step(mod0_birds, scope = formula(mod1_birds), perm.max = 200, na.action = na.exclude)
print(step_r2_birds)

save.image("emp_analysis_rda_2024_11_01a.RData")

# plot rda ----
# _ all hosts ----
df_metadata_sub = df_metadata_sub |> rownames_to_column() |> rename("sample_id" = "rowname")

# species
rda_scores_sp <- scores(rda_sp_all, display = "sites") # extract the site scores
rda_scores_sp_df <- as.data.frame(rda_scores_sp)
rda_scores_sp_df$sample_id <- rownames(rda_scores_sp_df) # add a study_id column
rda_fort_sp_all <- rda_scores_sp_df |> 
  left_join(df_metadata_sub, by = "sample_id") # join to metadata

# diet
rda_scores_diet <- scores(rda_diet_all, display = "sites")
rda_scores_diet_df <- as.data.frame(rda_scores_diet)
rda_scores_diet_df$sample_id <- rownames(rda_scores_diet_df)
rda_fort_diet_all <- rda_scores_diet_df |> 
  left_join(df_metadata_sub, by = "sample_id")

# sociality
rda_scores_social <- scores(rda_social_all, display = "sites")
rda_scores_social_df <- as.data.frame(rda_scores_social)
rda_scores_social_df$sample_id <- rownames(rda_scores_social_df)
rda_fort_social_all <- rda_scores_social_df |> 
  left_join(df_metadata_sub, by = "sample_id")

# plots
p_sp_all = ggplot(rda_fort_sp_all, aes(x = CAP1, y = CAP2, colour = host_species)) +
  geom_point() + theme(legend.position="none")
p_diet_all = ggplot(rda_fort_diet_all, aes(x = CAP1, y = CAP2, colour = basic_diet)) +
  geom_point()
p_social_all = ggplot(rda_fort_social_all, aes(x = CAP1, y = CAP2, colour = host_class)) +
  geom_point()+ theme(legend.position="none")

# _ mammals ----
df_metadata_sub_mammals = df_metadata_sub_mammals |> rownames_to_column() |> rename("sample_id" = "rowname")

# species
rda_scores_sp_mammals <- scores(rda_sp_mammals, display = "sites")
rda_scores_sp_mammals_df <- as.data.frame(rda_scores_sp_mammals)
rda_scores_sp_mammals_df$sample_id <- rownames(rda_scores_sp_mammals_df)
rda_fort_sp_mammals_all <- rda_scores_sp_mammals_df |> 
  left_join(df_metadata_sub_mammals, by = "sample_id")

# diet
rda_scores_diet_mammals <- scores(rda_diet_mammals, display = "sites")
rda_scores_diet_mammals_df <- as.data.frame(rda_scores_diet_mammals)
rda_scores_diet_mammals_df$sample_id <- rownames(rda_scores_diet_mammals_df)
rda_fort_diet_mammals_all <- rda_scores_diet_mammals_df |> 
  left_join(df_metadata_sub_mammals, by = "sample_id")

# sociality
rda_scores_social_mammals <- scores(rda_social_mammals, display = "sites")
rda_scores_social_mammals_df <- as.data.frame(rda_scores_social_mammals)
rda_scores_social_mammals_df$sample_id <- rownames(rda_scores_social_mammals_df)
rda_fort_social_mammals <- rda_scores_social_mammals_df |> 
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
rda_scores_sp_birds <- scores(rda_sp_birds, display = "sites") # extract the site scores
rda_scores_sp_df_birds <- as.data.frame(rda_scores_sp_birds)
rda_scores_sp_df_birds$sample_id <- rownames(rda_scores_sp_df_birds) # add a study_id column
rda_fort_sp_birds <- rda_scores_sp_df_birds |> 
  left_join(df_metadata_sub_birds, by = "sample_id") # join to metadata

# diet
rda_scores_diet_birds <- scores(rda_diet_birds, display = "sites")
rda_scores_diet_df_birds <- as.data.frame(rda_scores_diet_birds)
rda_scores_diet_df_birds$sample_id <- rownames(rda_scores_diet_df_birds)
rda_fort_diet_birds <- rda_scores_diet_df_birds |> 
  left_join(df_metadata_sub, by = "sample_id")

# sociality
rda_scores_social_birds <- scores(rda_social_birds, display = "sites")
rda_scores_social_df_birds <- as.data.frame(rda_scores_social_birds)
rda_scores_social_df_birds$sample_id <- rownames(rda_scores_social_df_birds)
rda_fort_social_birds <- rda_scores_social_df_birds |> 
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


# plots ----
par(mfrow = c(1,2))
aldex.plot(social_aldex_all, type = "MA", test = "welch", main = "MA plot")
aldex.plot(social_aldex_all, type = "MW", test = "welch", main = "effect plot")

# identify OTUs by name
taxonomy_table = p_tax |> as.data.frame() |> rownames_to_column()

# generate df of differentially abundant ASVs
social_diff_asvs = social_aldex_all |> as.data.frame() |> rownames_to_column() |> 
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

social_diff_asvs |> 
  ggplot(aes(x = effect, y = reorder(otu_scientific_name, (effect*-1)), colour = phylum)) +
  geom_point() +
  xlab("Effect size") +
  ylab("ASV") +
  theme(axis.text.y = element_text(face = "italic")) +
  ggtitle("Differentially abundant ASVs in Aves\nSocial and solitary\n(BH-corrected)") + 
  geom_vline(xintercept = 0, linetype = "dotted")

# _ birds ----
# generate df of differentially abundant ASVs
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
