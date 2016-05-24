if (require(testthat, quietly = TRUE)) {
  require(treatSens)
  test_check("treatSens")
} else {
  cat("package 'testthat' not available; cannot run unit tests\n")
}
