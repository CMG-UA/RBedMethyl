library(testthat)

if (requireNamespace("RBedMethyl", quietly = TRUE)) {
  library(RBedMethyl)
}

test_check("RBedMethyl")
