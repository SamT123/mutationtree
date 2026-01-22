# mutationtree

[![R-CMD-check](https://github.com/SamT123/mutationtree/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SamT123/mutationtree/actions/workflows/R-CMD-check.yaml)

Plot phylogenetic trees annotated with mutation annotations along branches.

## Installation

```r
devtools::install_github("SamT123/mutationtree")
```

## Relationship with `convergence`

`mutationtree` plots mutation annotated trees using the output of the [`convergence`](https://github.com/SamT123/convergence) package.

## Examples

### Basic Example

```r
library(tidyverse)
library(seqUtils)
library(convergence)
library(mutationtree)

# 1. Load and align sequences --------------------------------------------------
aligned_sequences <- "sequence.fasta" %>%
  seqUtils::fast_fasta() %>%
  seqUtils::mafft_align(seqUtils::alaska_232_2015_nts)

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
```

See `inst/example/basic_example.R` for a complete working example with H3N2 influenza sequences.

**Output:**

![Basic example output](inst/example/basic_output.png)

### Customized Example

Add annotations using four key arguments:

```r
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
```

See `inst/example/customized_example.R` for a complete working example.

**Output:**

![Customized example output](inst/example/customized_output.png)
