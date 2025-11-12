# Customized example: Advanced features for annotating mutation trees
# This example demonstrates lines, modified labels, and additional annotations

library(tidyverse)
library(seqUtils)
library(convergence)
library(mutationtree)

# Set paths for external tools (adjust these for your system)
cmaple_path <- "/Users/samturner/miniforge3/envs/treebuild/bin/"
usher_path <- "/Users/samturner/miniforge3/envs/intel_env/bin/"

# 1. Load and prepare data (same as basic example) -----------------------------
aligned_sequences <- "sequence.fasta" %>%
  seqUtils::fast_fasta() %>%
  seqUtils::mafft_align(seqUtils::alaska_232_2015_nts)

names(aligned_sequences) <- names(aligned_sequences) %>%
  stringr::str_remove_all(fixed(" ")) %>%
  stringr::str_split("[()]") %>%
  map_chr(~ paste0(.x[2:3], collapse = "_")) %>%
  trimws()

aligned_sequences <- aligned_sequences[
  !str_detect(names(aligned_sequences), "swine")
]

tree <- seqUtils::make_cmaple_tree(
  sequences = aligned_sequences,
  tree_path = "example.nwk",
  cmaple_path = cmaple_path,
  out_sequence = "A/Bayern/USAFSAM-16347/2025_H3N2",
  keep_files = c()
) %>%
  ape::ladderize()

tree_and_sequences <- convergence::makeTreeAndSequences(
  tree,
  tibble(
    Isolate_unique_identifier = names(aligned_sequences),
    dna_sequence = unname(aligned_sequences)
  )
)

tree_with_asr <- convergence::addASRusher(
  tree_and_sequences,
  aa_ref = seqUtils::alaska_232_2015_aas,
  nuc_ref = seqUtils::alaska_232_2015_nts,
  usher_path = usher_path
)

# 2. Find interesting nodes for annotation -------------------------------------
# Example: Find nodes with specific mutations
node_with_k189r <- tree_with_asr$tree_tibble %>%
  filter(map_lgl(aa_mutations_nonsyn, ~ "K189R" %in% .x)) %>%
  pull(node) %>%
  first()

# 3. Draw customized tree ------------------------------------------------------
mutationtree::draw_mutation_tree(
  tree_with_asr,
  file = "customized_output.png",
  width = 4,
  x_lim_expand = c(0.1, 0.5),
  # Add horizontal reference lines for clades
  lines = list(
    list(
      height = 21.5,
      text = "Clade of interest",
      line_color = "steelblue",
      text_color = "steelblue",
      text_cex = 0.8
    )
  ),

  # Highlight specific tip labels
  modified_tip_labels = list(
    list(
      old_label = "A/Bayern/USAFSAM-16347/2025_H3N2",
      new_label = "HIGHLIGHTED SEQUENCE: A/Bayern/USAFSAM-16347/2025_H3N2",
      text_color = "darkred",
      text_cex = 0.25
    )
  ),

  # Emphasize mutations at specific nodes (e.g., trunk)
  modified_node_labels = list(
    list(
      node = as.character(node_with_k189r),
      text_cex = 0.4,
      text_color = "purple"
    )
  ),

  # Add extra annotations with points
  additional_node_labels = list(
    list(
      node = node_with_k189r,
      text = "K189R emergence",
      text_color = "purple",
      point_color = "purple",
      point_cex = 1.0
    )
  )
)
