library(testthat)
library(devtools)

# Conditionally run rbundler tests, due to an issue with CRAN's policies:
# See https://github.com/opower/rbundler/issues/31
if(Sys.getenv('TEST_RBUNDLER') == TRUE) {
    message("Running Rbundler tests")
    test_check("rbundler")
} else {
    message("Not running Rbundler tests. Please set TEST_RBUNDLER=TRUE in your environment in order to run tests.")
}
