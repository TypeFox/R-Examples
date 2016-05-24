library(testthat)
library(devtools)
library(assertive.models)

with_envvar(
  c(LANG = "en_US"),
  test_check("assertive.models")
)
