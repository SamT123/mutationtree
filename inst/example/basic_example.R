# Basic example: Drawing a mutation-annotated phylogenetic tree
# This example uses H3N2 influenza HA sequences from GenBank

library(tidyverse)
library(seqUtils)
library(convergence)
library(mutationtree)

# Set paths for external tools (adjust these for your system)
cmaple_path <- "/Users/samturner/miniforge3/envs/treebuild/bin/"
usher_path <- "/Users/samturner/miniforge3/envs/intel_env/bin/"

# 1. Load and align sequences --------------------------------------------------
aligned_sequences <- "sequence.fasta" %>%
  seqUtils::fast_fasta() %>%
  seqUtils::mafft_align(seqUtils::alaska_232_2015_nts)

# Clean up sequence names
names(aligned_sequences) <- names(aligned_sequences) %>%
  stringr::str_remove_all(fixed(" ")) %>%
  stringr::str_split("[()]") %>%
  map_chr(~ paste0(.x[2:3], collapse = "_")) %>%
  trimws()

# Remove swine sequences
aligned_sequences <- aligned_sequences[
  !str_detect(names(aligned_sequences), "swine")
]

# 2. Build phylogenetic tree ---------------------------------------------------
tree <- seqUtils::make_cmaple_tree(
  sequences = aligned_sequences,
  tree_path = "example.nwk",
  cmaple_path = cmaple_path,
  out_sequence = "A/Bayern/USAFSAM-16347/2025_H3N2",
  keep_files = c()
) %>%
  ape::ladderize()

# 3. Create tree_and_sequences object ------------------------------------------
tree_and_sequences <- convergence::makeTreeAndSequences(
  tree,
  tibble(
    Isolate_unique_identifier = names(aligned_sequences),
    dna_sequence = unname(aligned_sequences)
  )
)

# 4. Add ancestral sequence reconstruction -------------------------------------
tree_with_asr <- convergence::addASRusher(
  tree_and_sequences,
  aa_ref = seqUtils::alaska_232_2015_aas,
  nuc_ref = seqUtils::alaska_232_2015_nts,
  usher_path = usher_path
)

# 5. Draw the mutation tree ----------------------------------------------------
mutationtree::draw_mutation_tree(
  tree_with_asr,
  file = "basic_output.png",
  width = 4,
  x_lim_expand = c(0.1, 0.5)
)
