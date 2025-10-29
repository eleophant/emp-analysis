########################################
#
# EMP analysis on subsampled data - host phylogeny construction
#
# Eleonore Lebeuf-Taylor
#
########################################

################ PACKAGES ################

library(ape)
library(rotl)
library(tidyverse)

################ DATA INPUT ################

load(file = "emp_data_for_host_phylogeny.RData")

################ SETUP ################

# 1) prepare species names (convert "sp." to genus level)
species_df <- tibble(
  original = unique(df_metadata_sub$host_species)
) %>%
  mutate(
    clean = str_remove(original, "^s__"),
    clean = str_replace_all(clean, "_", " "),
    clean_no_sp = str_remove(clean, " sp\\.")
  )

cat("Taxa to include in tree:\n")
print(species_df)

################ QUERY OTL ################

# 2) query Open Tree of Life
resolved_raw <- tnrs_match_names(species_df$clean_no_sp, context_name = "Animals")

# convert to tibble and join with original names
resolved <- as_tibble(resolved_raw) %>%
  mutate(search_string_lower = str_to_lower(search_string)) %>%
  left_join(
    species_df %>%
      mutate(clean_no_sp_lower = str_to_lower(clean_no_sp)),
    by = c("search_string_lower" = "clean_no_sp_lower")
  )

################ CHECKS ################

# 3) check what was matched and at what rank
cat("\nMatching summary:\n")
resolved %>%
  select(search_string, unique_name, approximate_match, score, original) %>%
  print(n = Inf)

# check all genus-level taxa
cat("\nAll genus-level matches:\n")
resolved %>%
  filter(str_detect(original, "_sp\\.")) %>%
  select(original, search_string, unique_name, ott_id, approximate_match) %>%
  print()

# check no more NAs in original column
cat("\nRows with NA in original column:\n")
resolved %>%
  filter(is.na(original)) %>%
  select(search_string, unique_name, ott_id) %>%
  print() # empty!

################ BUILD TREE ################

# 4) build tree
valid_ids <- resolved %>%
  filter(!is.na(ott_id)) %>%
  pull(ott_id)

cat(paste("Building tree with", length(valid_ids), "OTT IDs:\n"))
host_tree_complete <- tol_induced_subtree(ott_ids = valid_ids)

################ RENAME TIPS ################

# 5) rename tips
tip_to_original <- tibble(
  tree_tip_original = host_tree_complete$tip.label
) %>%
  mutate(
    tip_cleaned = tree_tip_original %>%
      str_remove("_ott[0-9]+") %>%
      str_replace_all("_", " ") %>%
      str_trim() %>%
      str_to_lower()
  ) %>%
  left_join(
    resolved %>%
      select(original, unique_name) %>%
      mutate(unique_lower = str_to_lower(str_trim(unique_name))),
    by = c("tip_cleaned" = "unique_lower")
  ) %>%
  mutate(renamed = coalesce(original, tree_tip_original))

host_tree_complete$tip.label <- tip_to_original$renamed

# check if Anser and Macropus are in tree
cat("\n=== Final Check ===\n")
cat("Anser_sp. in tree:", "s__Anser_sp." %in% host_tree_complete$tip.label, "\n")
cat("Macropus_sp. in tree:", "s__Macropus_sp." %in% host_tree_complete$tip.label, "\n")

# Anser is there, but Macropus is missing: add manually

################ ADD MACROPUS ################

cat("\nAdding Macropus_sp. manually:\n")

# list of potential relatives
macropus_relatives <- c(
  "Trichosurus",    # possums
  "Petaurus",       # gliders
  "Phascolarctos",  # koalas
  "Vombatus",       # wombats
  "Dasyurus",       # quolls
  "Sarcophilus",    # Tasmanian devil
  "Didelphis",      # opossums
  "Monodelphis",    # short-tailed opossums
  "Canis",          # dogs (other mammals)
  "Felis",          # cats
  "Vulpes",         # foxes
  "Mus",            # mice
  "Rattus"          # rats
)

# find a relative in the tree
relative_found <- NULL
for (genus in macropus_relatives) {
  matches <- grep(genus, host_tree_complete$tip.label, 
                  ignore.case = TRUE, value = TRUE)
  if (length(matches) > 0) {
    relative_found <- matches[1]
    cat(sprintf("  Found relative: %s\n", relative_found))
    break
  }
}

if (!is.null(relative_found)) {
  # add Macropus as sister to the relative
  tip_position <- which(host_tree_complete$tip.label == relative_found)
  
  host_tree_complete <- bind.tip(
    host_tree_complete,
    tip.label = "s__Macropus_sp.",
    where = tip_position,
    edge.length = 0.1,
    position = 0
  )
  
  cat(sprintf("  Added Macropus_sp. as sister to %s\n", relative_found))
  
}

################ FINAL CHECK ################

# check tree coverage
all_species <- unique(df_metadata_sub$host_species)
final_coverage <- sum(all_species %in% host_tree_complete$tip.label)
cat(sprintf("Coverage: %d/%d species (%.1f%%)\n",
            final_coverage, length(all_species),
            100 * final_coverage / length(all_species)))

# we have 100% coverage !

################ ADD BRANCHES ################

# 6) add branch lengths
host_phylo = host_tree_complete
host_phylo <- compute.brlen(host_phylo)