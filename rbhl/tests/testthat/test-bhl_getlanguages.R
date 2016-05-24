# tests for bhl_getlanguages fxn in rbhl
context("bhl_getlanguages")

test_that("bhl_getlanguages returns the correct class", {
  skip_on_cran()

  tt <- bhl_getlanguages()
  uu <- bhl_getlanguages('list')
  vv <- bhl_getlanguages('json')
  zz <- bhl_getlanguages('xml')

  expect_is(tt$data, "data.frame")

  expect_is(uu, "list")

  expect_is(vv, "character")
  expect_is(fromJSON(vv), "list")

  expect_is(zz, "character")
  expect_is(xmlParse(zz), "XMLInternalDocument")

  expect_equal(NCOL(tt$data), 2)
  expect_equal(length(uu), 3)
  expect_equal(length(uu$Status), 1)
  expect_equal(length(fromJSON(vv)), 3)
  expect_equal(length(zz), 1)
  expect_equal(length(xmlParse(zz)), 1)
})
