context("Pascal")

test_that("Pascal", {
  skip_on_cran()
  aux <- "../../inst/local-tests/local-Pascal.R"
  if(file.exists(aux)) source(aux, local=FALSE) else skip("skipped for R CMD Check")
  expect_identical(RHO[[4]], structure(c("0", "1/4", "1/2", "3/4", "1", "1/4", "0", "1/4",
                                         "1/2", "3/4", "1/2", "1/4", "0", "1/4", "1/2", "3/4", "1/2",
                                         "1/4", "0", "1/4", "1", "3/4", "1/2", "1/4", "0"), .Dim = c(5L,
                                                                                                     5L)))
})
