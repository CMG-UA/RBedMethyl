# RBedMethyl

Disk-backed access to nanoporetech modkit bedMethyl files for ONT-scale workflows.

## Installation

### Bioconductor (when available)

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("RBedMethyl")
```

### From source (local checkout)

```r
install.packages("devtools")
devtools::install(".")
```

## Usage

```r
library(RBedMethyl)

lines <- c(
  paste("chr1", 0, 1, "m", 0, "+", 0, 1, 0, 10, 0.5, 5, 5, 0, 0, 0, 0, 0, sep = "\t"),
  paste("chr1", 10, 11, "m", 0, "+", 10, 11, 0, 20, 0.25, 5, 15, 0, 0, 0, 0, 0, sep = "\t")
)

bed <- tempfile(fileext = ".bed")
writeLines(lines, bed)

h5 <- tempfile(fileext = ".h5")
bm <- readBedMethyl(bed, h5, mod = "m", fields = c("coverage", "pct", "mod_reads"))

# Per-site methylation fraction
beta(bm)

# Subset by region and summarize
regions <- GenomicRanges::GRanges(
  seqnames = "chr1",
  ranges = IRanges::IRanges(start = 1, end = 12)
)

bm_region <- subsetByRegion(bm, regions)
summarizeByRegion(bm, regions)

# Bracket subsetting with GRanges
bm[regions]
```
