context("rotation")

test_that("Rotation example", {
  skip_on_cran()
  aux <- "../../inst/local-tests/local-rotation.R"
  if(file.exists(aux)) source(aux, local=FALSE) else skip("skipped for R CMD Check")
  expect_identical(kanto2, as.bigq(1,9))
  expect_identical(kanto3, as.bigq(4,279))
})
