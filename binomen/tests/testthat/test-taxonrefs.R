context("taxonrefs")

test_that("taxonrefs basic functionality works", {
  a <- taxonref("genus", "Poa", 56, "http://scottchamberlain.info/")
  b <- taxonref("genus", "Quercus", 32343, "http://scottchamberlain.info/")
  aa <- taxonrefs(a, b)

  expect_is(aa, "taxonrefs")
  expect_is(aa[[1]], "taxonref")
  expect_is(aa[[2]], "taxonref")

  expect_equal(length(aa), 2)
  expect_null(names(aa))
})

test_that("taxonrefs check for class works", {
  expect_error(taxonrefs("4"), "One or more inputs was not of class taxonref")
  expect_error(taxonrefs(4), "One or more inputs was not of class taxonref")
  expect_error(taxonrefs(4, mtcars), "One or more inputs was not of class taxonref")
})
