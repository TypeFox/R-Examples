library(testthat)
library(threewords)

# Stolen from JennyB
if (identical(tolower(Sys.getenv("NOT_CRAN")), "true")) {
  test_check("threewords")
}