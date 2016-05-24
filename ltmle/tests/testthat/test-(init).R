context("Set seed")
test_that("Set seed", {
  set.seed(1)
  cat("seed set to 1\n")
  skip_on_cran() #keep seed constant on CRAN
  #skip("for testing as on cran - remove this to test everything!")
  set.seed(NULL)
  cat("seed set to NULL\n")
})
