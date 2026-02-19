testthat::test_that("readBedMethyl validates sorting", {
  lines <- c(
    paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr2", 0, 1, "m", 0, "-", 0, 1, 0, 8, 0.375, 3, 5, 0, 0, 0, 0, 0, sep = "\t")
  )

  tmp <- tempfile(fileext = ".bed")
  writeLines(lines, tmp)
  bm <- RBedMethyl::readBedMethyl(tmp, mod = "m", fields = c("coverage", "pct", "mod_reads"))
  testthat::expect_s4_class(bm, "RBedMethyl")
  testthat::expect_equal(as.numeric(RBedMethyl::beta(bm))[1], 0.5)

  unsorted <- tempfile(fileext = ".bed")
  writeLines(c(lines[2], lines[1]), unsorted)
  testthat::expect_error(
    RBedMethyl::readBedMethyl(unsorted, mod = "m", fields = c("coverage", "pct", "mod_reads")),
    "not sorted"
  )
})

testthat::test_that("summarizeByRegion and coercions work", {
  lines <- c(
    paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
  )

  tmp <- tempfile(fileext = ".bed")
  writeLines(lines, tmp)
  bm <- RBedMethyl::readBedMethyl(tmp, mod = "m", fields = c("coverage", "pct", "mod_reads"))

  regions <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 12)
  )

  summary <- RBedMethyl::summarizeByRegion(bm, regions)
  testthat::expect_equal(summary$coverage[1], 30)
  testthat::expect_equal(summary$mod_reads[1], 10)
  testthat::expect_equal(summary$beta[1], 10 / 30)

  rse <- as(bm, "RangedSummarizedExperiment")
  testthat::expect_s4_class(rse, "RangedSummarizedExperiment")

  if (requireNamespace("bsseq", quietly = TRUE) &&
      methods::hasMethod("coerce", c("RBedMethyl", "BSseq"))) {
    bs <- as(bm, "BSseq")
    testthat::expect_s4_class(bs, "BSseq")
  }
})

testthat::test_that("subsetByRegion supports GRanges and [ indexing", {
  lines <- c(
    paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr2", 0, 1, "m", 0, "-", 0, 1, 0, 8, 0.375, 3, 5, 0, 0, 0, 0, 0, sep = "\t")
  )

  tmp <- tempfile(fileext = ".bed")
  writeLines(lines, tmp)
  bm <- RBedMethyl::readBedMethyl(tmp, mod = "m", fields = c("coverage", "pct", "mod_reads"))

  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1, end = 5)
  )

  bm_gr <- RBedMethyl::subsetByRegion(bm, gr)
  testthat::expect_equal(length(bm_gr@index), 1L)
  testthat::expect_equal(as.numeric(RBedMethyl::beta(bm_gr))[1], 0.5)

  bm_bracket_gr <- bm[gr]
  testthat::expect_equal(length(bm_bracket_gr@index), 1L)

  bm_bracket_num <- bm[2]
  testthat::expect_equal(length(bm_bracket_num@index), 1L)
  testthat::expect_equal(as.numeric(RBedMethyl::beta(bm_bracket_num))[1], 0.25)

  bm_bracket_log <- bm[c(TRUE, FALSE, TRUE)]
  testthat::expect_equal(length(bm_bracket_log@index), 2L)
})

testthat::test_that("subsetByChromosomes filters by chromosome", {
  lines <- c(
    paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr2", 0, 1, "m", 0, "-", 0, 1, 0, 8, 0.375, 3, 5, 0, 0, 0, 0, 0, sep = "\t")
  )

  tmp <- tempfile(fileext = ".bed")
  writeLines(lines, tmp)
  bm <- RBedMethyl::readBedMethyl(tmp, mod = "m", fields = c("coverage", "pct", "mod_reads"))

  bm_chr1 <- RBedMethyl::subsetByChromosomes(bm, "chr1")
  testthat::expect_equal(length(bm_chr1@index), 2L)

  bm_chr2 <- RBedMethyl::subsetByChromosomes(bm, c("chr2"))
  testthat::expect_equal(length(bm_chr2@index), 1L)
})

testthat::test_that("example bedmethyl file can be read", {
  example_path <- system.file(
    "extdata",
    "example.bedmethyl",
    package = "RBedMethyl"
  )
  if (!file.exists(example_path)) {
    stop("example.bedmethyl is missing from the package extdata directory.")
  }

  bm <- RBedMethyl::readBedMethyl(example_path, mod = "m", fields = c("coverage", "pct", "mod_reads"))
  testthat::expect_s4_class(bm, "RBedMethyl")
  testthat::expect_true(length(bm@index) > 0L)
})

testthat::test_that("HDF5 file is reused on subsequent calls", {
  lines <- c(
    paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
  )

  tmp <- tempfile(fileext = ".bed")
  writeLines(lines, tmp)

  h5file <- tempfile(fileext = ".h5")

  # First call - creates the HDF5 file
  bm1 <- RBedMethyl::readBedMethyl(tmp, mod = "m", h5file = h5file,
                                    fields = c("coverage", "mod_reads"))
  testthat::expect_true(file.exists(h5file))
  mtime1 <- file.mtime(h5file)

  # Small delay to ensure mtime would differ if file were recreated
  Sys.sleep(0.1)

  # Second call - should reuse existing HDF5 file
  bm2 <- RBedMethyl::readBedMethyl(tmp, mod = "m", h5file = h5file,
                                    fields = c("coverage", "mod_reads"))
  mtime2 <- file.mtime(h5file)

  # File modification time should be unchanged (file was reused, not recreated)
  testthat::expect_equal(mtime1, mtime2)

  # Both objects should have identical data
  testthat::expect_equal(as.numeric(RBedMethyl::beta(bm1)), as.numeric(RBedMethyl::beta(bm2)))
  testthat::expect_equal(bm1@chrom_levels, bm2@chrom_levels)
  testthat::expect_equal(bm1@chr_index, bm2@chr_index)
})

testthat::test_that("HDF5 file is recreated when incomplete", {
  lines <- c(
    paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
  )

  tmp <- tempfile(fileext = ".bed")
  writeLines(lines, tmp)

  h5file <- tempfile(fileext = ".h5")

  # Create an incomplete/invalid HDF5 file
  rhdf5::h5createFile(h5file)
  rhdf5::h5write(1:10, h5file, "dummy_data")
  rhdf5::h5closeAll()
  mtime1 <- file.mtime(h5file)

  Sys.sleep(0.1)

  # Call should detect invalid file and recreate it
  bm <- RBedMethyl::readBedMethyl(tmp, mod = "m", h5file = h5file,
                                   fields = c("coverage", "mod_reads"))
  mtime2 <- file.mtime(h5file)

  # File should have been recreated (new mtime)
  testthat::expect_true(mtime2 > mtime1)
  testthat::expect_s4_class(bm, "RBedMethyl")
  testthat::expect_equal(length(bm@index), 2L)
})

testthat::test_that("HDF5 file is recreated when fields differ", {
  lines <- c(
    paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
  )

  tmp <- tempfile(fileext = ".bed")
  writeLines(lines, tmp)

  h5file <- tempfile(fileext = ".h5")

  # First call with coverage and mod_reads
  bm1 <- RBedMethyl::readBedMethyl(tmp, mod = "m", h5file = h5file,
                                    fields = c("coverage", "mod_reads"))
  testthat::expect_true(file.exists(h5file))
  mtime1 <- file.mtime(h5file)

  Sys.sleep(0.1)

  # Second call requesting additional field (pct) - should recreate
  bm2 <- RBedMethyl::readBedMethyl(tmp, mod = "m", h5file = h5file,
                                    fields = c("coverage", "mod_reads", "pct"))
  mtime2 <- file.mtime(h5file)

  # File should have been recreated to include pct field
  testthat::expect_true(mtime2 > mtime1)
  testthat::expect_true("pct" %in% names(bm2@assays))
})

testthat::test_that("show method prints a concise summary", {
  lines <- c(
    paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr2", 0, 1, "m", 0, "-", 0, 1, 0, 8, 0.375, 3, 5, 0, 0, 0, 0, 0, sep = "\t")
  )

  tmp <- tempfile(fileext = ".bed")
  writeLines(lines, tmp)
  bm <- RBedMethyl::readBedMethyl(tmp, mod = "m", fields = c("coverage", "mod_reads"))
  bm_chr1 <- RBedMethyl::subsetByChromosomes(bm, "chr1")

  out <- paste(capture.output(show(bm_chr1)), collapse = "\n")

  testthat::expect_match(out, "RBedMethyl object")
  testthat::expect_match(out, "Modification: m")
  testthat::expect_match(out, "Rows \\(active/total\\): 2 / 3")
  testthat::expect_match(out, "Chromosomes: 2")
  testthat::expect_match(out, "Assays: chrom, chromStart, chromEnd, strand, coverage, mod_reads")
})
