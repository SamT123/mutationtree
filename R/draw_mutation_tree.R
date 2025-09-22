#' Draw a mutation annotated tree
#'
#' @param tree_and_sequences_asr ...
#' @param file ...
#' @param width ...
#' @param height_per_sequence ...
#' @param lines ...
#' @param modified_tip_labels ...
#' @param modified_node_labels ...
#' @param additional_node_labels ...
#' @param tip_text_cex ...
#' @param tip_text_label_color ...
#' @param tip_text_nudge ...
#' @param tip_text_adj ...
#' @param node_text_cex ...
#' @param node_text_nudge ...
#' @param node_text_adj ...
#' @param tip_node_text_aa_mutation_color ...
#' @param tip_node_text_aa_mutation_special_colors ...
#' @param tip_node_text_syn_mutation_color ...
#' @param line_width ...
#' @param line_color ...
#' @param line_text_color ...
#' @param line_text_cex ...
#' @param line_text_adj ...
#' @param line_text_nudge ...
#' @param node_addtext_cex ...
#' @param node_addtext_color ...
#' @param node_addtext_nudge ...
#' @param node_addtext_adj ...
#' @param node_point_color ...
#' @param node_point_cex ...
#' @param edge_color ...
#' @param edge_width ...
#' @param x_lim_expand ...
#' @param y_lim_expand ...
#' @param node_nums ...
#'
#' @returns ...
#'
#' @importFrom dplyr filter
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
  node_text_nudge = c(-0.001, -0.001),
  node_text_adj = c(1, 1),

  tip_node_text_aa_mutation_color = "grey30",
  tip_node_text_aa_mutation_special_colors = list(
    red = c(145, 155, 156, 158, 159, 189, 193),
    orange = c(44L, 45L, 46L, 47L, 48L, 50L, 51L, 53L, 54L, 57L, 59L, 62L, 63L, 67L, 75L, 78L, 80L, 81L, 82L, 83L, 86L, 87L, 88L, 91L, 92L, 94L, 96L, 102L, 103L, 109L, 117L, 121L, 122L, 124L, 126L, 128L, 129L, 130L, 131L, 132L, 133L, 135L, 137L, 138L, 140L, 142L, 143L, 144L, 145L, 146L, 150L, 152L, 155L, 156L, 157L, 158L, 159L, 163L, 165L, 167L, 168L, 170L, 171L, 172L, 173L, 174L, 175L, 176L, 177L, 179L, 182L, 186L, 187L, 188L, 189L, 190L, 192L, 193L, 194L, 196L, 197L, 198L, 201L, 203L, 207L, 208L, 209L, 212L, 213L, 214L, 215L, 216L, 217L, 218L, 219L, 226L, 227L, 228L, 229L, 230L, 238L, 240L, 242L, 244L, 246L, 247L, 248L, 260L, 261L, 262L, 265L, 273L, 275L, 276L, 278L, 279L, 280L, 294L, 297L, 299L, 300L, 304L, 305L, 307L, 308L, 309L, 310L, 311L, 312L) # fmt: skip
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
      x0 = 0,
      x1 = 0 + 0.01 * x_lim[[2]] + 1,
      y0 = y_lim[[1]] + 0.05 * diff(y_lim),
      y1 = y_lim[[1]] + 0.05 * diff(y_lim),
      lwd = 0.7,
      col = "grey60"
    )
    text(
      x = 0.5,
      y = y_lim[[1]] + 0.05 * diff(y_lim),
      label = "1 mutation",
      cex = 0.5,
      adj = c(0.5, -0.5),
      col = "grey60"
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

      if (line_width > 0) {
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
    labbook::out.plot(
      code = draw(),
      fig_width = width,
      fig_height = height_per_sequence * ape::Ntip(tree_and_sequences_asr$tree)
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
