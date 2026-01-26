#' RBedMethyl class
#'
#' Disk-backed representation of nanoporetech modkit bedMethyl data from ONT sequencing.
#'
#' @slot assays A \code{SimpleList} of assay arrays.
#' @slot chrom_levels Character vector of chromosome names.
#' @slot strand_levels Character vector of strand levels.
#' @slot chr_index Matrix of chromosome row ranges (start/end).
#' @slot index Integer vector of active row indices.
#' @slot mod Modification code.
#'
#' @exportClass RBedMethyl
#' @importClassesFrom S4Vectors SimpleList
setClass(
  "RBedMethyl",
  slots = c(
    assays = "SimpleList",
    chrom_levels = "character",
    strand_levels = "character",
    chr_index = "matrix",
    index = "integer",
    mod = "character"
  )
)

setValidity("RBedMethyl", function(object) {
  required <- c("chrom", "chromStart", "chromEnd", "strand", "coverage")
  if (!all(required %in% names(object@assays))) {
    return("Missing required assay columns")
  }
  n <- length(object@assays$coverage)
  if (any(object@index < 1L | object@index > n)) {
    return("Index out of bounds")
  }
  TRUE
})
