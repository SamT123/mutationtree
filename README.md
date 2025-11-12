# mutationtree

[![R-CMD-check](https://github.com/SamT123/mutationtree/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SamT123/mutationtree/actions/workflows/R-CMD-check.yaml)

Plot phylogenetic trees annotated with amino acid and synonymous mutation annotations along branches.

## Installation

```r
devtools::install_github("SamT123/mutationtree")
```

## Relationship with convergence

`mutationtree` is a visualization layer built on top of [`convergence`](https://github.com/SamT123/convergence):

- **convergence** handles: tree building, ancestral sequence reconstruction, data structures
- **mutationtree** handles: plotting trees with mutation labels

## Examples

### Basic Example

Complete workflow from FASTA sequences to annotated tree:

```r
library(tidyverse)
library(seqUtils)
library(convergence)
library(mutationtree)

# 1. Load and align sequences
aligned_sequences <- "sequences.fasta" %>%
  seqUtils::fast_fasta() %>%
  seqUtils::mafft_align(seqUtils::alaska_232_2015_nts)

# 2. Build phylogenetic tree
tree <- seqUtils::make_cmaple_tree(
  sequences = aligned_sequences,
  tree_path = "tree.nwk",
  cmaple_path = "/path/to/cmaple/bin/",
  out_sequence = "outgroup_name"
) %>%
  ape::ladderize()

# 3. Create tree_and_sequences object
tree_and_sequences <- convergence::makeTreeAndSequences(
  tree,
  tibble(
    Isolate_unique_identifier = names(aligned_sequences),
    dna_sequence = unname(aligned_sequences)
  )
)

# 4. Add ancestral sequence reconstruction
tree_with_asr <- convergence::addASRusher(
  tree_and_sequences,
  aa_ref = seqUtils::alaska_232_2015_aas,
  nuc_ref = seqUtils::alaska_232_2015_nts,
  usher_path = "/path/to/usher/bin/"
)

# 5. Draw the tree
mutationtree::draw_mutation_tree(
  tree_with_asr,
  file = "output.png"
)
```

See `inst/example/basic_example.R` for a complete working example with H3N2 influenza sequences.

**Output:**

![Basic example output](inst/example/basic_output.png)

### Customized Example

Add annotations using four key arguments:

```r
mutationtree::draw_mutation_tree(
  tree_with_asr,
  file = "annotated.png",

  # Add horizontal reference lines for clades
  lines = list(
    list(
      height = 20,
      text = "Clade label",
      line_color = "steelblue",
      text_color = "steelblue"
    )
  ),

  # Highlight specific tips
  modified_tip_labels = list(
    list(
      old_label = "original_name",
      new_label = "Highlighted: original_name",
      text_color = "darkred",
      text_cex = 0.3
    )
  ),

  # Emphasize mutations at nodes
  modified_node_labels = list(
    list(
      node = "123",
      text_cex = 0.4,
      text_color = "purple"
    )
  ),

  # Add extra labels with points
  additional_node_labels = list(
    list(
      node = 123,
      text = "Important node",
      text_color = "purple",
      point_color = "purple"
    )
  )
)
```

See `inst/example/customized_example.R` for a complete working example.

**Output:**

![Customized example output](inst/example/customized_output.png)
