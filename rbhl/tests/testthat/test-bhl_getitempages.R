# tests for bhl_getitempages fxn in rbhl
context("bhl_getitempages")

test_that("bhl_getitempages returns the correct class", {
  skip_on_cran()

  tt <- bhl_getitempages('16800')
  uu <- bhl_getitempages('16800', as='xml')
  vv <- bhl_getitempages('16800', as='json')

  expect_is(tt$data, "data.frame")

  expect_is(bhl_getitempages('16800', as = "list")$Result, "list")

  expect_is(uu, "character")
  expect_is(xmlParse(uu), "XMLInternalDocument")
  expect_is(xpathApply(xmlParse(uu), "//Page"), "XMLNodeSet")

  expect_is(vv, "character")
  expect_is(jsonlite::fromJSON(vv, FALSE), "list")

  expect_equal(NCOL(tt$data), 13)
  expect_equal(length(uu), 1)
  expect_equal(length(xmlParse(uu)), 1)
  expect_equal(length(xpathApply(xmlParse(uu), "//Page")), 78)
  expect_equal(length(vv), 1)
  expect_equal(length(jsonlite::fromJSON(vv)), 3)
})
