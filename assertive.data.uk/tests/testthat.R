library(testthat)
library(devtools)
library(assertive.data.uk)

with_envvar(
  c(LANG = "en_US"),
  test_check("assertive.data.uk")
)
