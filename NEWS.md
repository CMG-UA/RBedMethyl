# RBedMethyl 0.99.0

- Initial Bioconductor submission with core functionality:
	- Disk-backed import of modkit bedMethyl files via `readBedMethyl()` with HDF5Array storage.
	- Field discovery with `bedMethylFields()`.
	- Per-site methylation fraction via `beta()`.
	- Assay-based filtering with `subsetBy()` and `filterByCoverage()`.
	- Genomic subsetting by interval and `GRanges` with `subsetByRegion()` and `[`.
	- Region-level summaries with `summarizeByRegion()`.
	- Coercion to `RangedSummarizedExperiment` and optional `BSseq`.
	- Vignette and runnable examples.
