library(testthat)
library(devtools)
library(assertive.matrices)

with_envvar(
  c(LANG = "en_US"),
  test_check("assertive.matrices")
)
