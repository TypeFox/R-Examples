library(testthat)
library(pingr)

if (Sys.getenv("NOT_CRAN") != "") {
  test_check("pingr")
}
