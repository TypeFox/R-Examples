# devtools::test(".", "calcDODN")
context("compDODN")

test_that("compDODN works as expected", {
  path <- system.file("extdata", package = "satellite")
  files <- list.files(path, pattern = glob2rx("LC8*.tif"), full.names = TRUE)
  sat <- satellite(files)

  t1 <- calcDODN(getSatDataLayer(sat, bcde = "B002n"))
  t2 <- calcDODN(getSatDataLayer(sat, bcde = "B004n"))
  expect_equal(t1, 8763)
  expect_equal(t2, 6677)
})
