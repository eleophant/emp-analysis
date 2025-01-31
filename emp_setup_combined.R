##***************************
## chapter 1: emp & sociability
##
## elt feat. kc
##
## 2023-11-24
##
##***************************

## startup packages ----
# packages for converting biom to df
library(tidyverse)
library(biomformat)
library(conflicted)

# packages for analysis
library(microbiome)
library(vegan)
library(viridis)
theme_set(theme_bw())

## biom to df conversion ----

emp = read_biom("emp_cr_silva_16S_123.subset_2k.biom")

## manipulations to use in R
# convert the biom (sparse) into a matrix (dense = zeroes are added in)
df_otu_emp = 
  emp |> 
  biom_data() |> 
  as.matrix() |> 
  t() |> # transpose; necessary to follow standard statistical conventions
  as_tibble(rownames = "sample_id") # keep sample id
# Warning message: In asMethod(object) :
#sparse->dense coercion: allocating vector of size 1.9 GiB
df_otu_emp
write_tsv(df_otu_emp, "df_otu_emp.tsv")

df_sample_emp_full = read_tsv("emp_qiime_mapping_qc_filtered.tsv") |>
  rename(sample_id = '#SampleID')
df_sample_emp
write_tsv(df_sample_emp_full, "df_sample_emp_full.tsv")

df_taxa_emp = emp |> 
  observation_metadata() |> 
  as_tibble(rownames = "taxa_id")
df_taxa_emp
write_tsv(df_taxa_emp, "df_taxa_emp.tsv")

# link sample metadata to otu data (possible because they're both tibbles, not matrices)
df_vert_gut_dense = df_sample_emp |> 
  right_join(df_otu_emp) |> # keep order of rows from 2nd df (df_otu_emp)
  filter(sample_scientific_name == "gut metagenome" & host_phylum == "p__Chordata" & host_scientific_name != "Homo sapiens")

vc_sums = df_vert_gut_dense |> 
  select(77:126806) |> 
  colSums() # yields sum for each column
vc_sums

sum(vc_sums == 0) # tells you how many taxa are completely absent from subset

# exclude taxa that are 0 in all samples (absent from subset)
df_vert_gut = df_vert_gut_dense |> 
  select(1:76, names(which(vc_sums != 0)))

write_tsv(df_vert_gut, "df_vert_gut.tsv")

df_taxa_gut = df_taxa_emp |> 
  filter(taxa_id %in% names(which(vc_sums != 0))) |> 
  #remove chloroplast and mitochondrial DNA
  filter(taxonomy3 != "D_2__Chloroplast") |>
  filter(taxonomy5 != "D_4__Mitochondria")
write_tsv(df_taxa_gut_full, "df_taxa_gut_full.tsv")

# check that there are no humans in vertebrate gut metagenome df
which(df_vert_gut$host_scientific_name == "Homo sapiens") #result: integer(0)

# from OTU taxonomy table, make a list of all OTUs that are bacteria
vc_bacteria = df_taxa_gut |> 
  filter(taxonomy1 == "D_0__Bacteria") |> 
  pull(taxa_id)

# remove non-bacterial OTUs from df_vert_gut if desired
df_vert_gut |> 
  select(1:76, all_of(vc_bacteria))
write_tsv(df_vert_gut, "df_vert_gut.tsv")

# end of startup

## initial data input ----
df_otu_emp = read_tsv("df_otu_emp_full.tsv") #result of converting biom object into a dense (with 0s) matrix
df_sample_emp = read_tsv("df_sample_emp_full.tsv") #metadata (mapping file)
df_taxa_emp = read_tsv("df_taxa_emp_full.tsv") #taxonomic table of all OTUs
df_vert_gut = read_tsv("df_vert_gut_full.tsv") #metadata + OTU table of all vertebrate gut metagenome samples

## fix metadata names ----
df_vert_gut = df_vert_gut |> 
  # identify NAs in host_common_name
  #filter(is.na(host_common_name)) |> 
  #pull(host_species) #|> 
  # fix missing host_common_name
  mutate(host_common_name = replace(host_common_name, host_species == "s__Carollia_sowelli", "Sowell's short-tailed bat")) |> 
  mutate(host_common_name = replace(host_common_name, host_species == "s__Thraupis_palmarum", "palm tanager")) |> 
  mutate(host_common_name = replace(host_common_name, host_species == "s__Turdus_olivater", "black-hooded thrush")) |>
  mutate(host_common_name = replace(host_common_name, host_species == "s__Ursus_americanus", "black bear")) |>
  # fix missing host_class for Iguana iguana
  mutate(host_class = replace(host_class, host_class == "c__", "c__Reptilia")) |> 
  #fix missing host_order for Rusa unicolor
  mutate(host_order = replace(host_order, host_scientific_name == "Rusa unicolor", "o__Artidactyla")) #|>
#check:
#filter(host_order == "o__") |> 
#pull(host_scientific_name)

## exclusions ----
df_vert_gut = df_vert_gut |> 
  filter(!is.na(host_common_name)) |> 
  # remove colobine primates (captivity)
  filter(host_scientific_name != "Pygathrix nemaeus") |>
  # remove hosts not identified at species
  filter(host_species != "s__") |> 
  # keep only lower intestine samples
  filter(!grepl(c("Crop"), Description)) |> 
  filter(!grepl(c("Gizzard"), Description)) |> 
  filter(!grepl(c("Upper"), Description)) |> 
  #remove duplicate samples
  filter(!duplicated(host_subject_id))

## checks
which(df_vert_gut$host_common_name == "Francois langur")
which(df_vert_gut$host_species == "Trachypithecus francoisi")
which(df_vert_gut$host_common_name == "Silvered-leaf langur")
which(df_vert_gut$host_species == "Trachypithecus cristatus")
which(df_vert_gut$host_common_name_provided == "Northern Douc")
#integer(0), removal worked

write_tsv(df_vert_gut, "df_vert_gut.tsv")


## merge social data ----

# read in df with sociality data
host_sociality = read_csv("host_species_social_traits.csv")
host_sociality = host_sociality |> 
  select(host_scientific_name, basic_sociality)
#host_sociality = unique(host_sociality)

# remove any duplicate species
host_sociality = host_sociality |> 
  arrange(host_scientific_name, basic_sociality) |> 
  filter(duplicated(host_scientific_name) == FALSE)

# map the basic_sociality column to species found in metadata
metadata = df_vert_gut |> 
  left_join(host_sociality, by = "host_scientific_name")

colnames(metadata) #basic_sociality is there!

## write tsv files ----
# create bacterial community file
community = df_vert_gut |> 
  select(1, 77:56402)
write_tsv(community, "community.tsv")

# create metadata file
metadata = df_vert_gut |> 
  select(1:76)
write_tsv(metadata, "metadata.tsv")

# build phyloseq ----
community = read_tsv("../full_emp_files/community.tsv")
metadata = read_tsv("metadata.tsv")
otu_taxonomy = read_tsv("df_taxa_gut_full.tsv")
tree_emp = read_tree("silva_123.97_otus.tre", errorIfNULL = TRUE)

# _ otu table ----
# make an otu_table with correct row names
community_pivoted = community |> 
  pivot_longer(-sample_id) |> 
  pivot_wider(names_from = sample_id)
# coerce to a df and make the first columns = row names (tibbles don't accept row names)
community_pivoted = as.data.frame(community_pivoted)
row.names(community_pivoted) = community_pivoted[,1]
community_pivoted = community_pivoted[,-1]

# now try to create a phyloseq otu table from the transposed community df
p_otu = otu_table(community_pivoted, taxa_are_rows = TRUE)
sample_names(p_otu) # it works!


# _metadata table ----
metadata = metadata |> 
  mutate(basic_sociality = replace(basic_sociality, host_scientific_name == "Ursus americanus americanus", "solitary")) |> 
  # add basic_diet column
  mutate(basic_diet = case_when(
    diet=="carnivorous" ~ "carnivorous",
    diet=="blood" ~ "carnivorous",
    diet=="insectivorous" ~ "carnivorous",
    diet=="piscivorous" ~ "carnivorous",
    diet=="myrmecophagus" ~ "carnivorous",
    diet=="myrmecophagus-frugivorous" ~ "omnivorous",
    diet=="frugivore-insectivore" ~ "omnivorous",
    diet=="frugivorous-insectivorous" ~ "omnivorous",
    diet=="herbivorous" ~ "herbivorous",
    diet=="granivorous" ~ "herbivorous",
    diet=="frugivorous" ~ "herbivorous",
    diet=="nectarivorous" ~ "herbivorous",
    diet=="omnivorous" ~ "omnivorous",
    diet=="omnivorous" ~ "omnivorous",
  ))

p_metadata = sample_data(metadata)

# attempt at making metadata with sample_id as row names
# currently, sample_names(p_metadata_social) returns sa1, sa2, ...
# try assign first column as row names
row.names(p_metadata) = p_metadata$sample_id
#sample_names(p_metadata_social) # matches sample_names(p_otu)!


# _tax_table ----

# first, rename the columns of df_taxa_gut to taxonomic levels
head(otu_taxonomy)
otu_taxonomy = otu_taxonomy |> 
  dplyr::rename(domain = taxonomy1, phylum = taxonomy2, class = taxonomy3, order = taxonomy4, family = taxonomy5, genus = taxonomy6, species = taxonomy7)
# convert to df (can't set row names on a tibble)
otu_taxonomy = as.data.frame(otu_taxonomy)
row.names(otu_taxonomy) = otu_taxonomy$taxa_id
otu_taxonomy = otu_taxonomy[,-1, drop = FALSE]

# now create a phyloseq tax_table
p_tax = tax_table(as.matrix(otu_taxonomy)) #works!


## create phyloseq ----
physeq = phyloseq(p_otu, p_metadata, p_tax) #works!

# merge tree
ps = merge_phyloseq(physeq, tree_emp)
# everything works! we now have a phyloseq object ready to be used for our analysis

# assign sample_id as row names in metadata; currently, sample_names(p_metadata_social) returns sa1, sa2, so need to assign first column as row names
p_metadata = sample_data(metadata)
row.names(p_metadata) = p_metadata$sample_id
head(sample_names(p_metadata)) 
#check if it matches p_otu sample names
head(sample_names(p_otu)) #it works!

## phyloseq checks
sample_names(ps) #looks good!
rank_names(ps) #"domain"  "phylum"  "class"   "order"   "family"  "genus"   "species"
sample_variables(ps) #all correct!
head(otu_table(ps))
head(tax_table(ps))
head(taxa_names(ps))
phy_tree(tree_emp)
## everything checks out!

## updating phyloseq metadata ----
# change metadata to df with basic_diet column
p_metadata = sample_data(metadata)
row.names(p_metadata) = p_metadata$sample_id
sample_data(ps) = p_metadata
colnames(sample_data(ps)) #basic_diet is there!

## subsetting phyloseq object ----
ps_mammals = subset_samples(ps, host_class == "c__Mammalia")
ps_birds = subset_samples(ps, host_class == "c__Aves")
ps_phylum = tax_glom(ps, taxrank = "phylum", NArm = FALSE)

