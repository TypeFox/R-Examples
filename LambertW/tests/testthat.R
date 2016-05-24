library(testthat)

cat("Starting tests for 'LambertW' package in 1 second ...")
Sys.sleep(1)

set.seed(10231)

test_check("LambertW")
