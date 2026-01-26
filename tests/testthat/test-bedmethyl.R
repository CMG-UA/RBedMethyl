testthat::test_that("readBedMethyl validates sorting", {
  lines <- c(
    paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t"),
    paste("chr2", 0, 1, "m", 0, "-", 0, 1, 0, 8, 0.375, 3, 5, 0, 0, 0, 0, 0, sep = "\t")
  )

  tmp <- tempfile(fileext = ".bed")
  writeLines(lines, tmp)
  h5 <- tempfile(fileext = ".h5")

  bm <- RBedMethyl::readBedMethyl(tmp, h5, mod = "m", fields = c("coverage", "pct", "mod_reads"))
  testthat::expect_s4_class(bm, "RBedMethyl")
  testthat::expect_equal(as.numeric(RBedMethyl::beta(bm))[1], 0.5)

  unsorted <- tempfile(fileext = ".bed")
  writeLines(c(lines[2], lines[1]), unsorted)
  testthat::expect_error(
    RBedMethyl::readBedMethyl(unsorted, tempfile(fileext = ".h5"), mod = "m", fields = c("coverage", "pct", "mod_reads")),
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
  h5 <- tempfile(fileext = ".h5")
  bm <- RBedMethyl::readBedMethyl(tmp, h5, mod = "m", fields = c("coverage", "pct", "mod_reads"))

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

  if (requireNamespace("bsseq", quietly = TRUE)) {
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
  h5 <- tempfile(fileext = ".h5")
  bm <- RBedMethyl::readBedMethyl(tmp, h5, mod = "m", fields = c("coverage", "pct", "mod_reads"))

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

testthat::test_that("example bedmethyl file can be read", {
  example_path <- file.path(
    system.file(package = "RBedMethyl"),
    "data",
    "example.bedmethyl"
  )
  if (!file.exists(example_path)) {
    stop("example.bedmethyl is missing from the package data directory.")
  }

  h5 <- tempfile(fileext = ".h5")
  bm <- RBedMethyl::readBedMethyl(example_path, h5, mod = "m", fields = c("coverage", "pct", "mod_reads"))
  testthat::expect_s4_class(bm, "RBedMethyl")
  testthat::expect_true(length(bm@index) > 0L)
})
