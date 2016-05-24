# Unsetting R_TESTS is a workaround for an R bug, see
# http://github.com/hadley/testthat/issues/144
# Should be removed once the issue is resolved
Sys.setenv("R_TESTS" = "")
library(testthat)
library(mkin)

test_check("mkin")
