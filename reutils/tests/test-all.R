library(testthat)
opts <- options(reutils.verbose.queries = FALSE, reutils.email = "test@mail.com")
test_check("reutils")
options(opts)