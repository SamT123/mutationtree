#' Draw a mutation-annotated phylogenetic tree
#'
#' Creates a visualization of a phylogenetic tree with non-synonymous (amino acid)
#' and synonymous mutations labeled on branches. Supports extensive customization
#' of labels, colors, and annotations.
#'
#' @param tree_and_sequences_asr A tree_and_sequences object with ancestral sequence
#'   reconstruction, created by \code{convergence::addASRusher()}.
#' @param file Output file path. If NULL, displays interactively. Supported formats:
#'   PNG (600 dpi) or PDF.
#' @param width Figure width in inches. Default: 6.
#' @param height_per_sequence Height per tip in inches. Total height is calculated
#'   as \code{height_per_sequence * n_tips}. Default: 0.03.
#' @param lines List of lists defining horizontal reference lines. Each element:
#'   \code{list(height, text, line_color, text_color, line_width, text_cex, text_adj, text_nudge)}.
#' @param modified_tip_labels List of lists to customize tip labels. Each element:
#'   \code{list(old_label, new_label, text_color, text_cex)}.
#' @param modified_node_labels List of lists to customize node mutation labels. Each element:
#'   \code{list(node, text_color, text_cex)}.
#' @param additional_node_labels List of lists to add extra annotations. Each element:
#'   \code{list(node, text, text_color, text_cex, text_adj, text_nudge, point_color, point_cex)}.
#' @param tip_text_cex Tip label text size. Default: 0.2.
#' @param tip_text_label_color Tip label text color. Default: "grey70".
#' @param tip_text_nudge Tip label position offset as \code{c(x, y)} fraction of plot range. Default: \code{c(0.002, 0)}.
#' @param tip_text_adj Tip label text justification. Default: \code{c(0, 0.5)}.
#' @param node_text_cex Node mutation label text size. Default: 0.2.
#' @param node_text_nudge Node label position offset as \code{c(x, y)} fraction of plot range. Default: \code{c(-0.002, -0.015)}.
#' @param node_text_adj Node label text justification. Default: \code{c(1, 1)}.
#' @param tip_node_text_aa_mutation_color Default color for amino acid mutations. Default: "grey30".
#' @param tip_node_text_aa_mutation_special_colors Named list mapping colors to amino acid positions
#'   to highlight, e.g. \code{list(red = c(145, 155), blue = c(189))}. Default: H3 epitope positions.
#' @param tip_node_text_syn_mutation_color Color for synonymous mutation counts. Default: "grey70".
#' @param line_width Width of horizontal reference lines. Default: 0.3.
#' @param line_color Color of horizontal reference lines. Default: "grey50".
#' @param line_text_color Color of reference line text. Default: "black".
#' @param line_text_cex Reference line text size. Default: 1.
#' @param line_text_adj Reference line text justification. Default: \code{c(1, 1)}.
#' @param line_text_nudge Reference line text position offset. Default: \code{c(-0.01, -0.003)}.
#' @param node_addtext_cex Additional node label text size. Default: \code{node_text_cex * 4}.
#' @param node_addtext_color Additional node label text color. Default: "black".
#' @param node_addtext_nudge Additional node label position offset. Default: \code{c(-0.004, -0.006)}.
#' @param node_addtext_adj Additional node label text justification. Default: \code{c(1, 1)}.
#' @param node_point_color Additional node point color. Default: "black".
#' @param node_point_cex Additional node point size. Default: 0.8.
#' @param edge_color Color of tree edges. Default: "grey80".
#' @param edge_width Width of tree edges. If NULL, uses ape default. Default: NULL.
#' @param x_lim_expand Plot x-axis expansion as \code{c(left, right)} fractions. Default: \code{c(0.02, 0.2)}.
#' @param y_lim_expand Plot y-axis expansion as \code{c(bottom, top)} fractions. Default: \code{c(0.02, 0.02)}.
#' @param node_nums Logical; show node numbers for debugging. Default: FALSE.
#'
#' @returns If \code{file} is specified, returns the file path. If \code{file} is NULL,
#'   returns a list with components: \code{draw_fn} (function to execute the plot),
#'   \code{width} (plot width), and \code{height} (plot height).
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' draw_mutation_tree(tree_with_asr)
#'
#' # Save to file
#' draw_mutation_tree(tree_with_asr, file = "tree.png")
#'
#' # With annotations
#' draw_mutation_tree(
#'   tree_with_asr,
#'   lines = list(list(height = 20, text = "Clade A")),
#'   modified_tip_labels = list(list(old_label = "seq1", text_color = "red"))
#' )
#' }
#'
#' @importFrom dplyr filter
#' @importFrom graphics segments
#' @export
draw_mutation_tree = function(
  tree_and_sequences_asr,
  file = NULL,
  width = 6,
  height_per_sequence = 0.03,
  lines = list(),
  modified_tip_labels = list(),
  modified_node_labels = list(),
  additional_node_labels = list(),

  # default parameters
  tip_text_cex = 0.2, # text_cex
  tip_text_label_color = "grey70",
  tip_text_nudge = c(0.002, 0),
  tip_text_adj = c(0, 0.5),

  node_text_cex = 0.2,
  node_text_nudge = c(-0.002, -0.015),
  node_text_adj = c(1, 1),

  tip_node_text_aa_mutation_color = "grey30",
  tip_node_text_aa_mutation_special_colors = list(
    red = seqUtils::h3_epitope_positions[["koel"]],
    orange = seqUtils::h3_epitope_positions[["wolf"]]
  ),
  tip_node_text_syn_mutation_color = "grey70",

  line_width = 0.3,
  line_color = "grey50",
  line_text_color = "black",
  line_text_cex = 1,
  line_text_adj = c(1, 1),
  line_text_nudge = c(-0.01, -0.003),

  node_addtext_cex = node_text_cex * 4,
  node_addtext_color = "black",
  node_addtext_nudge = c(-0.004, -0.006),
  node_addtext_adj = c(1, 1),

  node_point_color = "black",
  node_point_cex = 0.8,

  edge_color = "grey80",
  edge_width = NULL, # NULL -> use ape default

  x_lim_expand = c(0.02, 0.2),
  y_lim_expand = c(0.02, 0.02),
  node_nums = F
) {
  draw_tree_inner = function(
    tree_and_sequences_asr,
    lines,
    modified_tip_labels,
    modified_node_labels,
    additional_node_labels,

    # default parameters
    tip_text_cex,
    tip_text_label_color,
    tip_text_nudge,
    tip_text_adj,

    node_text_cex,
    node_text_nudge,
    node_text_adj,

    tip_node_text_aa_mutation_color,
    tip_node_text_aa_mutation_special_colors,
    tip_node_text_syn_mutation_color,

    line_width,
    line_color,
    line_text_color,
    line_text_cex,
    line_text_adj,
    line_text_nudge,

    node_addtext_cex,
    node_addtext_color,
    node_addtext_nudge,
    node_addtext_adj,

    node_point_color,
    node_point_cex,

    edge_color,
    edge_width,

    x_lim_expand,
    y_lim_expand,
    node_nums
  ) {
    # set up plot ---
    par(mai = rep(0, 4), xaxs = "i", yaxs = "i")

    x_lim = range(ape::node.depth.edgelength(tree_and_sequences_asr$tree))
    x_lim = x_lim + diff(x_lim) * x_lim_expand * c(-1, 1)

    y_lim = range(ape::node.height(tree_and_sequences_asr$tree))
    y_lim = y_lim + diff(y_lim) * y_lim_expand * c(-1, 1)

    # plot bare tree ---
    ape::plot.phylo(
      tree_and_sequences_asr$tree,
      show.tip.label = F,
      edge.color = edge_color,
      edge.width = edge_width,
      no.margin = T,
      x.lim = x_lim,
      y.lim = y_lim
    )

    segments(
      x0 = 0 + 0.01 * x_lim[[2]],
      x1 = 0 + 0.01 * x_lim[[2]] + 0.001,
      y0 = y_lim[[1]] + 0.1 * diff(y_lim),
      y1 = y_lim[[1]] + 0.1 * diff(y_lim),
      lwd = 0.7,
      col = "grey40"
    )
    text(
      x = 0 + 0.01 * x_lim[[2]] + 0.001 / 2,
      y = y_lim[[1]] + 0.1 * diff(y_lim),
      label = "1/1000 s/s",
      cex = 0.3,
      adj = c(0.5, -0.8),
      col = "grey40"
    )

    # node mutation labels ---

    set_if_present = function(l, nm, default) {
      if (nm %in% names(l)) {
        return(l[[nm]])
      } else {
        return(default)
      }
    }

    internal_nodes = filter(
      tree_and_sequences_asr$tree_tibble,
      node %in% parent
    )$node

    modified_node_labels = setNames(
      modified_node_labels,
      purrr::map_chr(modified_node_labels, "node")
    )

    stopifnot(!any(duplicated(names(modified_node_labels)))) # do not modify the same node label twice

    for (node in internal_nodes) {
      if (as.character(node) %in% names(modified_node_labels)) {
        mn = modified_node_labels[[as.character(node)]]
        # fmt: skip
        this_node_color = set_if_present(mn, "text_color", tip_node_text_aa_mutation_color)
        this_node_cex = set_if_present(mn, "text_cex", node_text_cex)

        if (
          length(tip_node_text_aa_mutation_special_colors) > 0 &
            "text_color" %in% mn
        ) {
          warning(
            "`modified_node_label` text_color overridden by `tip_node_text_aa_mutation_special_colors`"
          )
        }
      } else {
        this_node_color = tip_node_text_aa_mutation_color
        this_node_cex = node_text_cex
      }

      labels = get_labels_for_one_node(
        tree_and_sequences_asr$tree_tibble$aa_mutations_nonsyn[[node]],
        length(tree_and_sequences_asr$tree_tibble$nt_mutations_syn[[node]]),
        aa_position_colors = tip_node_text_aa_mutation_special_colors,
        default_aa_color = this_node_color,
        syn_color = tip_node_text_syn_mutation_color
      )

      text(
        x = ape::node.depth.edgelength(tree_and_sequences_asr$tree)[[node]] +
          node_text_nudge[[1]] * diff(x_lim),
        y = ape::node.height(tree_and_sequences_asr$tree)[[node]] +
          node_text_nudge[[2]] * diff(y_lim),
        label = labels$labels,
        col = labels$colors,
        adj = node_text_adj,
        cex = this_node_cex,
        family = font_family
      )

      if (node_nums) {
        text(
          x = ape::node.depth.edgelength(tree_and_sequences_asr$tree)[[node]],
          y = ape::node.height(tree_and_sequences_asr$tree)[[node]],
          label = tree_and_sequences_asr$tree_tibble$node[[node]],
          col = "black",
          adj = c(-0.2, 0.5),
          cex = node_text_cex,
          family = font_family
        )
      }
    }

    # tip name & mutation labels ---
    tips = filter(
      tree_and_sequences_asr$tree_tibble,
      !node %in% parent
    )$node

    modified_tip_labels = setNames(
      modified_tip_labels,
      purrr::map_chr(modified_tip_labels, "old_label")
    )

    stopifnot(!any(duplicated(names(modified_tip_labels)))) # do not modify the same node label twice

    for (node in tips) {
      this_tip_label = tree_and_sequences_asr$tree_tibble$label[[node]]

      if (this_tip_label %in% names(modified_tip_labels)) {
        mt = modified_tip_labels[[this_tip_label]]
        this_tip_color = set_if_present(mt, "text_color", tip_text_label_color)
        this_tip_cex = set_if_present(mt, "text_cex", tip_text_cex)
        this_tip_label_mod = set_if_present(mt, "new_label", this_tip_label)
      } else {
        this_tip_color = tip_text_label_color
        this_tip_cex = tip_text_cex
        this_tip_label_mod = this_tip_label
      }

      labels = get_labels_for_one_tip(
        this_tip_label_mod,
        tree_and_sequences_asr$tree_tibble$aa_mutations_nonsyn[[node]],
        length(tree_and_sequences_asr$tree_tibble$nt_mutations_syn[[node]]),
        aa_position_colors = tip_node_text_aa_mutation_special_colors,
        default_aa_color = tip_node_text_aa_mutation_color,
        syn_color = tip_node_text_syn_mutation_color,
        tip_label_color = this_tip_color
      )

      text(
        x = ape::node.depth.edgelength(tree_and_sequences_asr$tree)[[node]] +
          tip_text_nudge[[1]] * diff(x_lim),
        y = ape::node.height(tree_and_sequences_asr$tree)[[node]] +
          tip_text_nudge[[2]] * diff(y_lim),
        label = labels$labels,
        col = labels$colors,
        adj = tip_text_adj,
        cex = this_tip_cex,
        family = font_family
      )
    }

    # horizontal lines ---

    ## this is the easiest way to set the defaults
    draw_line = function(
      height,
      line_width = get("line_width", envir = parent.frame()),
      line_color = get("line_color", envir = parent.frame()),
      text = NA,
      text_color = line_text_color,
      text_cex = line_text_cex,
      text_adj = line_text_adj,
      text_nudge = line_text_nudge
    ) {
      if (line_width > 0) {
        abline(
          h = height,
          lw = line_width,
          col = line_color
        )
      }

      if (!is.na(text)) {
        text(
          label = text,
          y = height + text_nudge[[2]] * diff(y_lim),
          x = x_lim[[2]] + text_nudge[[1]] * diff(x_lim),
          col = text_color,
          cex = text_cex,
          adj = text_adj
        )
      }
    }

    for (line in lines) {
      do.call(what = draw_line, args = line)
    }

    # additional node labels ---

    ## this is the easiest way to set the defaults
    draw_additional_node_label = function(
      node,
      text = NA,
      text_cex = node_addtext_cex,
      text_color = node_addtext_color,
      text_nudge = node_addtext_nudge,
      text_adj = node_addtext_adj,
      point_color = node_point_color,
      point_cex = node_point_cex
    ) {
      x = ape::node.depth.edgelength(tree_and_sequences_asr$tree)[[node]]
      y = ape::node.height(tree_and_sequences_asr$tree)[[node]]

      if (point_cex > 0) {
        points(
          x = x,
          y = y,
          pch = 16,
          cex = point_cex,
          col = point_color
        )
      }

      if (!is.na(text)) {
        text(
          label = text,
          x = x + text_nudge[[1]] * diff(x_lim),
          y = y + text_nudge[[2]] * diff(y_lim),
          col = text_color,
          cex = text_cex,
          adj = text_adj
        )
      }
    }

    for (additional_node in additional_node_labels) {
      do.call(
        what = draw_additional_node_label,
        args = additional_node
      )
    }
  }

  draw = function() {
    draw_tree_inner(
      tree_and_sequences_asr = tree_and_sequences_asr,
      lines = lines,
      modified_tip_labels = modified_tip_labels,
      modified_node_labels = modified_node_labels,
      additional_node_labels = additional_node_labels,
      # default parameter
      tip_text_cex = tip_text_cex,
      tip_text_label_color = tip_text_label_color,
      tip_text_nudge = tip_text_nudge,
      tip_text_adj = tip_text_adj,

      node_text_cex = node_text_cex,
      node_text_nudge = node_text_nudge,
      node_text_adj = node_text_adj,

      tip_node_text_aa_mutation_color = tip_node_text_aa_mutation_color,
      tip_node_text_aa_mutation_special_colors = tip_node_text_aa_mutation_special_colors,
      tip_node_text_syn_mutation_color = tip_node_text_syn_mutation_color,

      line_width = line_width,
      line_color = line_color,
      line_text_color = line_text_color,
      line_text_cex = line_text_cex,
      line_text_adj = line_text_adj,
      line_text_nudge = line_text_nudge,

      node_addtext_cex = node_addtext_cex,
      node_addtext_color = node_addtext_color,
      node_addtext_nudge = node_addtext_nudge,
      node_addtext_adj = node_addtext_adj,

      node_point_color = node_point_color,
      node_point_cex = node_point_cex,

      edge_color = edge_color,
      edge_width = edge_width,

      x_lim_expand = x_lim_expand,
      y_lim_expand = y_lim_expand,
      node_num = node_nums
    )
  }

  if (!"Droid Sans Mono" %in% sysfonts::font_families()) {
    message(
      "Droid Sans Mono font not found, so using default mono font. Add it with sysfonts::font_add & showtext::showtext_auto."
    )
    font_family = "mono"
  } else {
    font_family = "Droid Sans Mono"
  }

  # file
  if (is.null(file)) {
    return(
      list(
        draw_fn = draw,
        width = width,
        height = height_per_sequence * ape::Ntip(tree_and_sequences_asr$tree)
      )
    )
  } else {
    format = fs::path_ext(file)
    stopifnot(format %in% c("png", "pdf"))

    if (format == "pdf") {
      pdf(
        file = file,
        width = width,
        height = height_per_sequence * ape::Ntip(tree_and_sequences_asr$tree)
      )
    } else if (format == "png") {
      png(
        filename = file,
        width = width,
        height = height_per_sequence * ape::Ntip(tree_and_sequences_asr$tree),
        units = "in",
        res = 600
      )
    }

    draw()

    graphics.off()
    return(file)
  }
}

## making labels ----------------------------------------
get_mutation_labels = function(
  aa_mutations_nonsyn,
  num_mutations_syn,
  aa_position_colors,
  direction,
  default_aa_color,
  syn_color
) {
  if (length(aa_position_colors) > 0) {
    aa_position_colors = purrr::imap(
      aa_position_colors,
      ~ setNames(rep(.y, length(.x)), as.character(.x))
    ) %>%
      purrr::reduce(.f = c)
  } else {
    aa_position_colors = NULL
  }

  add_label_or_color = function(
    existing_labels,
    new_label,
    direction,
    pad = T
  ) {
    npad = ifelse(
      length(existing_labels) > 0,
      nchar(existing_labels[[1]]) + 1,
      0
    )

    if (pad) {
      padding = paste0(rep(" ", npad), collapse = "")
    } else {
      padding = ""
    }

    if (direction == "left") {
      existing_labels = c(
        paste0(new_label, padding),
        existing_labels
      )
    } else {
      existing_labels = c(
        paste0(padding, new_label),
        existing_labels
      )
    }

    existing_labels
  }

  labels = c()
  colors = c()

  if (direction == "left") {
    aa_mutations_nonsyn = rev(aa_mutations_nonsyn)
  }

  # add aa mutations
  for (aa_mutation in aa_mutations_nonsyn) {
    labels = add_label_or_color(labels, aa_mutation, direction)

    at = stringr::str_sub(aa_mutation, 2, -2)

    if (at %in% names(aa_position_colors)) {
      colors = add_label_or_color(
        colors,
        aa_position_colors[[at]],
        direction,
        pad = F
      )
    } else {
      colors = add_label_or_color(
        colors,
        default_aa_color,
        direction,
        pad = F
      )
    }
  }

  # add syn nts
  if (num_mutations_syn > 0) {
    labels = add_label_or_color(
      labels,
      paste0("+", num_mutations_syn),
      direction
    )
    colors = add_label_or_color(
      colors,
      syn_color,
      direction,
      pad = F
    )
  }

  list(labels = labels, colors = colors)
}

get_labels_for_one_node = function(
  aa_mutations_nonsyn,
  num_mutations_syn,
  aa_position_colors,
  default_aa_color,
  syn_color
) {
  get_mutation_labels(
    aa_mutations_nonsyn,
    num_mutations_syn,
    aa_position_colors,
    direction = "left",
    default_aa_color = default_aa_color,
    syn_color = syn_color
  )
}

get_labels_for_one_tip = function(
  tip_label,
  aa_mutations_nonsyn,
  num_mutations_syn,
  aa_position_colors,
  default_aa_color,
  syn_color,
  tip_label_color
) {
  labels = get_mutation_labels(
    aa_mutations_nonsyn,
    num_mutations_syn,
    aa_position_colors,
    direction = "right",
    default_aa_color = default_aa_color,
    syn_color = syn_color
  )

  labels$labels = paste0(
    paste0(rep(" ", nchar(tip_label) + 1), collapse = ""),
    labels$labels
  )

  labels$labels = c(tip_label, labels$labels)
  labels$colors = c(tip_label_color, labels$colors)

  labels
}
