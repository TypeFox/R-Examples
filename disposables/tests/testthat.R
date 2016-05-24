
if (Sys.getenv("NOT_CRAN") != "") {
  library(testthat)
  library(disposables)
  test_check("disposables")
}
