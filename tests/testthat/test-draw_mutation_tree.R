test_that("draw_mutation_tree returns correct output", {
  skip_if_not_installed("seqUtils")
  skip_if_not_installed("convergence")

  # Get example data
  example_fasta <- system.file(
    "example",
    "sequence.fasta",
    package = "mutationtree"
  )
  skip_if(example_fasta == "", "Example FASTA not found")

  # Skip if external tools are not available
  cmaple_path <- "/Users/samturner/miniforge3/envs/treebuild/bin/"
  usher_path <- "/Users/samturner/miniforge3/envs/intel_env/bin/"

  skip_if_not(file.exists(file.path(cmaple_path, "cmaple")), "cmaple not found")
  skip_if_not(file.exists(file.path(usher_path, "usher")), "usher not found")

  # Run basic example workflow to get tree_with_asr
  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)

  aligned_sequences <- example_fasta %>%
    seqUtils::fast_fasta() %>%
    seqUtils::mafft_align(seqUtils::alaska_232_2015_nts)

  names(aligned_sequences) <- names(aligned_sequences) %>%
    stringr::str_remove_all(stringr::fixed(" ")) %>%
    stringr::str_split("[()]") %>%
    purrr::map_chr(~ paste0(.x[2:3], collapse = "_")) %>%
    trimws()

  aligned_sequences <- aligned_sequences[
    !stringr::str_detect(names(aligned_sequences), "swine")
  ]

  tree <- seqUtils::make_cmaple_tree(
    sequences = aligned_sequences,
    tree_path = "test_tree.nwk",
    cmaple_path = cmaple_path,
    out_sequence = "A/Bayern/USAFSAM-16347/2025_H3N2",
    keep_files = c()
  ) %>%
    ape::ladderize()

  tree_and_sequences <- convergence::makeTreeAndSequences(
    tree,
    tibble::tibble(
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

  # 5. Test file output
  result <- draw_mutation_tree(
    tree_with_asr,
    file = "test_output.png"
  )

  expect_equal(result, "test_output.png")
  expect_true(file.exists("test_output.png"))

  # 6. Test without file output
  result_no_file <- draw_mutation_tree(
    tree_with_asr,
    file = NULL
  )

  expect_type(result_no_file, "list")
  expect_true("draw_fn" %in% names(result_no_file))
  expect_true("width" %in% names(result_no_file))
  expect_true("height" %in% names(result_no_file))

  # Clean up
  setwd(old_wd)
  if (file.exists(file.path(temp_dir, "test_output.png"))) {
    unlink(file.path(temp_dir, "test_output.png"))
  }
})
