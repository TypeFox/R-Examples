library(testthat)
library(CHNOSZ)

# as fix for https://github.com/hadley/testthat/issues/129, https://github.com/hadley/testthat/issues/144
# (since we use makeCluster() etc via palply)
Sys.setenv("R_TESTS" = "") 

test_check("CHNOSZ")
