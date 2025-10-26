library(ape)
library(rotl)
library(tidyverse)

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

# Verify no more NAs in original column
cat("\nRows with NA in original column:\n")
resolved %>%
  filter(is.na(original)) %>%
  select(search_string, unique_name, ott_id) %>%
  print()

# 4) build tree
valid_ids <- resolved %>%
  filter(!is.na(ott_id)) %>%
  pull(ott_id)

cat(paste("Building tree with", length(valid_ids), "OTT IDs...\n"))
host_tree_complete <- tol_induced_subtree(ott_ids = valid_ids)

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
)

# check tree coverage
all_species <- unique(df_metadata_sub$host_species)
final_coverage <- sum(all_species %in% host_tree_complete$tip.label)
cat(sprintf("Coverage: %d/%d species (%.1f%%)\n",
            final_coverage, length(all_species),
            100 * final_coverage / length(all_species)))

# we have 100% coverage !

host_phylo = host_tree_complete
