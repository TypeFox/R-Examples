# tests for bhl_getpanames fxn in rbhl
context("bhl_getsubjecttitles")

test_that("bhl_getsubjecttitles returns the correct class", {
	skip_on_cran()

	tt <- bhl_getsubjecttitles('diptera')
	uu <- bhl_getsubjecttitles('diptera', 'xml')
	vv <- bhl_getsubjecttitles('diptera', 'json')

  expect_is(tt$data, "data.frame")

  expect_is(uu, "character")
  expect_is(xmlParse(uu), "XMLInternalDocument")

  expect_is(vv, "character")
  expect_is(jsonlite::fromJSON(vv), "list")

  expect_equal(length(uu), 1)
  expect_equal(length(xmlParse(uu)), 1)
  expect_equal(length(vv), 1)
  expect_equal(length(jsonlite::fromJSON(vv)), 3)
})
