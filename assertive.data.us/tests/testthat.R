library(testthat)
library(devtools)
library(assertive.data.us)

with_envvar(
  c(LANG = "en_US"),
  test_check("assertive.data.us")
)
