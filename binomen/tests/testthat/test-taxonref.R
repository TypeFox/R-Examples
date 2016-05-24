context("taxonref")

test_that("taxonref basic functionality works", {
  aa <- taxonref("genus", "Poa", 56, "http://scottchamberlain.info/")

  expect_is(aa, "taxonref")
  expect_is(aa$rank, "character")
  expect_is(aa$name, "character")
  expect_is(aa$id, "character")
  expect_is(aa$uri, "character")

  expect_equal(aa$rank, "genus")
  expect_equal(aa$id, "56")
})

test_that("taxonref with no input", {
  aa <- taxonref()

  expect_is(aa, "taxonref")
  expect_is(aa$rank, "character")
  expect_is(aa$name, "character")
  expect_is(aa$id, "character")
  expect_is(aa$uri, "character")

  expect_equal(aa$rank, "none")
  expect_equal(aa$id, "none")
})
