context("Testing distribution family type\n")

location.scale <- c("normal", "cauchy", "t", "unif")
scale.only <- c("gamma", "f", "exp", "chisq")

test_that("location/scale/non-negativity are correct", {
  for (dd in location.scale) {
    expect_identical(get_distname_family(dd),
                     list(location = TRUE,
                          scale = TRUE,
                          is.non.negative = FALSE))
  }
  
  for (dd in scale.only) {
    expect_identical(get_distname_family(dd),
                     list(location = FALSE,
                          scale = TRUE,
                          is.non.negative = TRUE))
  }
})