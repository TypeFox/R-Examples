# tests for bhl_booksearch fxn in rbhl
context("bhl_booksearch")

test_that("bhl_booksearch returns the correct class", {
	skip_on_cran()

	tt <- bhl_booksearch('evolution', year=2000)
	uu <- bhl_booksearch('evolution', year=2000, as='list')
	vv <- bhl_booksearch('evolution', year=2000, as="xml")

  expect_is(tt$data, "data.frame")

	expect_is(uu, "list")

	expect_is(vv, "character")
	expect_is(xmlParse(vv), "XMLInternalDocument")
	expect_is(xpathApply(xmlParse(vv), '//Result'), "XMLNodeSet")
})
