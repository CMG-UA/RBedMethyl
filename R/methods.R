if (!isGeneric("beta")) {
  setGeneric("beta", function(a, b) standardGeneric("beta"))
}

#' Show an RBedMethyl summary
#'
#' Print a concise summary of an \code{RBedMethyl} object.
setMethod("show", "RBedMethyl", function(object) {
  total_rows <- length(object@assays$chrom)
  active_rows <- length(object@index)
  assays <- names(object@assays)

  n_chr <- length(object@chrom_levels)
  n_preview <- min(5L, n_chr)
  chrom_preview <- if (n_preview > 0L) {
    paste(object@chrom_levels[seq_len(n_preview)], collapse = ", ")
  } else {
    "<none>"
  }
  if (n_chr > n_preview) {
    chrom_preview <- paste0(chrom_preview, ", ...")
  }

  cat("RBedMethyl object\n")
  cat("  Modification:", paste(object@mod, collapse = ", "), "\n")
  cat("  Rows (active/total):", active_rows, "/", total_rows, "\n")
  cat("  Chromosomes:", n_chr, "(", chrom_preview, ")\n")
  cat("  Assays:", paste(assays, collapse = ", "), "\n")
})

#' Per-site methylation fraction
#'
#' Compute per-site methylation fraction for an \code{RBedMethyl} object.
#' Requires the \code{mod_reads} assay to be loaded.
#'
#' @param a An \code{RBedMethyl} object.
#' @param b Unused, kept for \code{base::beta} compatibility.
#'
#' @return Numeric vector of per-site methylation fractions.
#' @export
#' @importFrom methods is validObject
#' @aliases beta,RBedMethyl-method
setMethod("beta", signature(a = "RBedMethyl", b = "missing"), function(a, b) {
  if ("pct" %in% names(a@assays)) {
    return(a@assays$pct[a@index])
  }
  if (!"mod_reads" %in% names(a@assays)) {
    stop("mod_reads assay not available; load it via fields.")
  }
  a@assays$mod_reads[a@index] / a@assays$coverage[a@index]
})


#' Subset by assay predicate
#'
#' Subset an \code{RBedMethyl} object using a predicate over an assay.
#'
#' @param x An \code{RBedMethyl} object.
#' @param column Assay name to filter on (must be loaded).
#' @param FUN Predicate function returning a logical vector.
#'
#' @return A filtered \code{RBedMethyl} object.
#' @export
#' @aliases subsetBy,RBedMethyl-method
#' @examples
#' lines <- c(
#'   paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
#'   paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
#' )
#' tmp <- tempfile(fileext = ".bed")
#' writeLines(lines, tmp)
#' bm <- readBedMethyl(tmp, mod = "m", fields = c("coverage", "pct", "mod_reads"))
#' bm2 <- subsetBy(bm, "coverage", function(v) v >= 15)
#' length(RBedMethyl::beta(bm2))
setGeneric("subsetBy", function(x, column, FUN) standardGeneric("subsetBy"))
#' @export
setMethod("subsetBy", "RBedMethyl", function(x, column, FUN) {
  if (!column %in% names(x@assays)) {
    stop("Assay not available: ", column)
  }
  vals <- x@assays[[column]][x@index]
  keep <- which(FUN(vals))
  x@index <- x@index[keep]
  validObject(x)
  x
})

#' Filter by coverage
#'
#' Filter an \code{RBedMethyl} object by minimum coverage.
#'
#' @param x An \code{RBedMethyl} object.
#' @param min_cov Minimum coverage threshold.
#'
#' @return A filtered \code{RBedMethyl} object.
#' @export
#' @aliases filterCoverage,RBedMethyl-method
#' @examples
#' lines <- c(
#'   paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
#'   paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
#' )
#' tmp <- tempfile(fileext = ".bed")
#' writeLines(lines, tmp)
#' bm <- readBedMethyl(tmp, mod = "m", fields = c("coverage", "pct", "mod_reads"))
#' bm2 <- filterCoverage(bm, min_cov = 15)
#' length(RBedMethyl::beta(bm2))
setGeneric("filterCoverage", function(x, min_cov) standardGeneric("filterCoverage"))
#' @export
setMethod("filterCoverage", "RBedMethyl", function(x, min_cov = 5L) {
  subsetBy(x, "coverage", function(v) v >= min_cov)
})

#' Subset by chromosomes
#'
#' Subset an \code{RBedMethyl} object by one or more chromosomes.
#'
#' @param x An \code{RBedMethyl} object.
#' @param chr Character vector of chromosome names.
#'
#' @return A filtered \code{RBedMethyl} object.
#' @export
#' @aliases subsetByChromosomes,RBedMethyl-method
#' @examples
#' lines <- c(
#'   paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
#'   paste("chr2", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
#' )
#' tmp <- tempfile(fileext = ".bed")
#' writeLines(lines, tmp)
#' bm <- readBedMethyl(tmp, mod = "m", fields = c("coverage", "pct", "mod_reads"))
#' bm2 <- subsetByChromosomes(bm, c("chr1"))
#' length(RBedMethyl::beta(bm2))
setGeneric("subsetByChromosomes", function(x, chr) standardGeneric("subsetByChromosomes"))
#' @export
setMethod("subsetByChromosomes", "RBedMethyl", function(x, chr) {
  if (!is.character(chr) || length(chr) == 0L) {
    stop("chr must be a non-empty character vector.")
  }
  chr <- unique(chr)
  missing_chr <- setdiff(chr, rownames(x@chr_index))
  if (length(missing_chr) > 0L) {
    stop("Requested chromosome(s) not present in object: ", paste(missing_chr, collapse = ", "))
  }

  starts <- vapply(chr, function(cn) x@chr_index[cn, 1], integer(1L))
  ends <- vapply(chr, function(cn) x@chr_index[cn, 2], integer(1L))
  keep <- unlist(Map(seq.int, starts, ends), use.names = FALSE)

  x@index <- intersect(x@index, keep)
  validObject(x)
  x
})


#' Subset by region
#'
#' Subset an \code{RBedMethyl} object by genomic interval.
#'
#' @param x An \code{RBedMethyl} object.
#' @param chr Chromosome name.
#' @param start Region start (0-based, half-open).
#' @param end Region end.
#'
#' @return A filtered \code{RBedMethyl} object.
#' @export
#' @aliases subsetByRegion,RBedMethyl-method
#' @aliases subsetByRegion,RBedMethyl,character,numeric,numeric-method
#' @examples
#' lines <- c(
#'   paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
#'   paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
#' )
#' tmp <- tempfile(fileext = ".bed")
#' writeLines(lines, tmp)
#' bm <- readBedMethyl(tmp, mod = "m", fields = c("coverage", "pct", "mod_reads"))
#' bm2 <- subsetByRegion(bm, "chr1", 0, 5)
#' length(RBedMethyl::beta(bm2))
setGeneric(
  "subsetByRegion",
  function(x, chr, start, end) standardGeneric("subsetByRegion")
)

#' @export
setMethod(
  "subsetByRegion",
  signature(x = "RBedMethyl", chr = "character", start = "numeric", end = "numeric"),
  function(x, chr, start, end) {
    if (length(chr) != 1L) {
      stop("chr must be a single chromosome.")
    }
    if (!chr %in% rownames(x@chr_index)) {
      stop("Requested chromosome not present in object.")
    }
    rows <- seq.int(x@chr_index[chr, 1], x@chr_index[chr, 2])
    keep <- rows[
      x@assays$chromEnd[rows] > start &
        x@assays$chromStart[rows] < end
    ]
    x@index <- intersect(x@index, keep)
    validObject(x)
    x
  }
)

#' Subset by GRanges
#'
#' Subset an \code{RBedMethyl} object by overlaps with a \code{GRanges}.
#'
#' @param x An \code{RBedMethyl} object.
#' @param chr A \code{GRanges} object of regions.
#' @param start Unused (for signature compatibility).
#' @param end Unused (for signature compatibility).
#'
#' @return A filtered \code{RBedMethyl} object.
#' @export
#' @aliases subsetByRegion,RBedMethyl,GRanges-method
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
setMethod(
  "subsetByRegion",
  signature(x = "RBedMethyl", chr = "GRanges", start = "missing", end = "missing"),
  function(x, chr, start, end) {
    regions <- chr
    if (length(x@index) == 0L) {
      return(x)
    }

    idx <- x@index
    chrom <- x@chrom_levels[x@assays$chrom[idx]]
    start_pos <- x@assays$chromStart[idx]
    end_pos <- x@assays$chromEnd[idx]
    strand <- x@strand_levels[x@assays$strand[idx]]

    gr <- GenomicRanges::GRanges(
      seqnames = chrom,
      ranges = IRanges::IRanges(start = start_pos + 1L, end = end_pos),
      strand = strand
    )

    hits <- GenomicRanges::findOverlaps(regions, gr, ignore.strand = TRUE)
    rows <- unique(S4Vectors::subjectHits(hits))
    x@index <- idx[rows]
    validObject(x)
    x
  }
)

#' Subset rows
#'
#' Subset an \code{RBedMethyl} object by integer, logical, or \code{GRanges} index.
#'
#' @param x An \code{RBedMethyl} object.
#' @param i Integer, logical, or \code{GRanges} index.
#' @param j Unused.
#' @param ... Unused.
#' @param drop Unused.
#'
#' @return A filtered \code{RBedMethyl} object.
#' @export
#' @aliases [,RBedMethyl-method
#' @aliases [,RBedMethyl,missing,missing,missing-method
#' @aliases [,RBedMethyl,ANY,missing,missing-method
setMethod("[", signature(x = "RBedMethyl", i = "missing", j = "missing", drop = "missing"),
  function(x, i, j, ..., drop) {
    x
  }
)

#' @export
setMethod("[", signature(x = "RBedMethyl", i = "ANY", j = "missing", drop = "missing"),
  function(x, i, j, ..., drop) {
    if (is(i, "GRanges")) {
      return(subsetByRegion(x, i))
    }

    if (is.logical(i)) {
      if (length(i) != length(x@index)) {
        stop("Logical index length must match number of active rows.")
      }
      x@index <- x@index[i]
    } else if (is.numeric(i) || is.integer(i)) {
      x@index <- x@index[i]
    } else {
      stop("Unsupported index type for RBedMethyl: ", class(i)[1])
    }

    validObject(x)
    x
  }
)

#' Summarize by regions
#'
#' Summarize methylation by a set of regions.
#'
#' @param x An \code{RBedMethyl} object.
#' @param regions A \code{GRanges} of regions.
#'
#' @return A \code{DataFrame} with coverage, mod_reads, beta, and n_sites.
#' @export
#' @aliases summarizeByRegion,RBedMethyl-method
#' @examples
#' lines <- c(
#'   paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
#'   paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
#' )
#' tmp <- tempfile(fileext = ".bed")
#' writeLines(lines, tmp)
#' bm <- readBedMethyl(tmp, mod = "m", fields = c("coverage", "pct", "mod_reads"))
#' regions <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges::IRanges(start = 1, end = 12)
#' )
#' summarizeByRegion(bm, regions)
setGeneric(
  "summarizeByRegion",
  function(x, regions) standardGeneric("summarizeByRegion")
)

#' @export
setMethod("summarizeByRegion", "RBedMethyl", function(x, regions) {
  if (!is(regions, "GRanges")) {
    stop("regions must be a GRanges object.")
  }
  if (!"mod_reads" %in% names(x@assays)) {
    stop("mod_reads assay not available; load it via fields.")
  }
  if (!"coverage" %in% names(x@assays)) {
    stop("coverage assay not available; load it via fields.")
  }

  idx <- x@index
  chrom <- x@chrom_levels[x@assays$chrom[idx]]
  start <- x@assays$chromStart[idx]
  end <- x@assays$chromEnd[idx]
  strand <- x@strand_levels[x@assays$strand[idx]]

  gr <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(start = start + 1L, end = end),
    strand = strand
  )

  hits <- GenomicRanges::findOverlaps(regions, gr, ignore.strand = TRUE)
  by_region <- split(S4Vectors::subjectHits(hits), S4Vectors::queryHits(hits))

  coverage <- x@assays$coverage[idx]
  mod_reads <- x@assays$mod_reads[idx]

  region_ids <- seq_along(regions)
  cov_sum <- numeric(length(regions))
  mod_sum <- numeric(length(regions))
  n_sites <- integer(length(regions))

  for (i in region_ids) {
    rows <- by_region[[as.character(i)]]
    if (is.null(rows) || length(rows) == 0L) {
      cov_sum[i] <- 0
      mod_sum[i] <- 0
      n_sites[i] <- 0L
    } else {
      mat <- DelayedArray::cbind(coverage[rows], mod_reads[rows])
      sums <- DelayedMatrixStats::colSums2(mat)
      cov_sum[i] <- sums[1]
      mod_sum[i] <- sums[2]
      n_sites[i] <- length(rows)
    }
  }

  beta <- ifelse(cov_sum > 0, mod_sum / cov_sum, NA_real_)
  S4Vectors::DataFrame(
    coverage = cov_sum,
    mod_reads = mod_sum,
    beta = beta,
    n_sites = n_sites,
    row.names = names(regions)
  )
})

setAs("RBedMethyl", "RangedSummarizedExperiment", function(from) {
  idx <- from@index
  chrom <- from@chrom_levels[from@assays$chrom[idx]]
  start <- from@assays$chromStart[idx]
  end <- from@assays$chromEnd[idx]
  strand <- from@strand_levels[from@assays$strand[idx]]

  gr <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(start = start + 1L, end = end),
    strand = strand
  )

  assays <- S4Vectors::SimpleList()
  if ("coverage" %in% names(from@assays)) {
    assays$coverage <- DelayedArray::cbind(from@assays$coverage[idx])
  }
  if ("mod_reads" %in% names(from@assays)) {
    assays$mod_reads <- DelayedArray::cbind(from@assays$mod_reads[idx])
  }
  if ("pct" %in% names(from@assays)) {
    assays$pct <- DelayedArray::cbind(from@assays$pct[idx])
  }

  SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    rowRanges = gr
  )
})

if (requireNamespace("bsseq", quietly = TRUE)) {
  setAs("RBedMethyl", "BSseq", function(from) {
    if (!all(c("coverage", "mod_reads") %in% names(from@assays))) {
      stop("coverage and mod_reads assays are required for BSseq coercion.")
    }
    idx <- from@index
    chrom <- from@chrom_levels[from@assays$chrom[idx]]
    start <- from@assays$chromStart[idx]
    end <- from@assays$chromEnd[idx]
    strand <- from@strand_levels[from@assays$strand[idx]]

    gr <- GenomicRanges::GRanges(
      seqnames = chrom,
      ranges = IRanges::IRanges(start = start + 1L, end = end),
      strand = strand
    )

    M <- DelayedArray::cbind(from@assays$mod_reads[idx])
    Cov <- DelayedArray::cbind(from@assays$coverage[idx])

    bsseq::BSseq(M = M, Cov = Cov, gr = gr)
  })
}
