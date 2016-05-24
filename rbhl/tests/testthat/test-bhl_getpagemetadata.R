# tests for bhl_getpagemetadata fxn in rbhl
context("bhl_getpagemetadata")

test_that("bhl_getpagemetadata returns the correct class", {
	skip_on_cran()

	tt <- bhl_getpagemetadata(page=1328690, ocr=TRUE, as='list')
	uu <- bhl_getpagemetadata(page=1328690, ocr=TRUE, as='xml')
	vv <- bhl_getpagemetadata(page=1328690, ocr=TRUE, as='json')

  expect_is(tt, "list")
  expect_is(tt$Result, "list")
  expect_null(tt$ErrorMessage)

  expect_is(uu, "character")

  expect_is(vv, "character")


  expect_equal(length(tt), 3)
  expect_equal(length(tt$Status), 1)
  expect_equal(length(uu), 1)
  expect_equal(length(vv), 1)
})
