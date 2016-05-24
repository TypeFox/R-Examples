# tests for bhl_getitemmetadata fxn in rbhl
context("bhl_getitemmetadata")

test_that("bhl_getitemmetadata returns the correct class", {
	skip_on_cran()

	tt <- bhl_getitemmetadata('16800', TRUE)
	uu <- bhl_getitemmetadata('16800', TRUE, as='xml')
	vv <- bhl_getitemmetadata('16800', TRUE, as='json')

  expect_is(tt$data, "data.frame")

  expect_is(uu, "character")
  expect_is(xmlParse(uu), "XMLInternalDocument")

  expect_is(vv, "character")
  expect_is(fromJSON(vv), "list")

  expect_equal(NCOL(tt$data), 31)
  expect_equal(length(uu), 1)
  expect_equal(length(vv), 1)
  expect_equal(length(fromJSON(vv)), 3)
})
