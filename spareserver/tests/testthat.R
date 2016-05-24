library(testthat)
library(spareserver)

if (Sys.getenv("NOT_CRAN") != "") {
  test_check("spareserver")
}
