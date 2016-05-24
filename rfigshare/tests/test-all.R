if (packageVersion("testthat") >= "0.7.1.99") {
    library(testthat)
  test_check("rfigshare")
} else {
  library(testthat)
  test_package("rfigshare")
}
