library(testthat)
library(withr)
library(assertive.sets)

with_options(
  c(useFancyQuotes = FALSE),
  test_check("assertive.sets")
)
