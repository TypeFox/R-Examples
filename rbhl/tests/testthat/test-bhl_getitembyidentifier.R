# tests for bhl_getitembyidentifier fxn in rbhl
context("bhl_getitembyidentifier")

test_that("bhl_getitembyidentifier returns the correct class", {
	skip_on_cran()

	tt <- bhl_getitembyidentifier(type='ia', value='animalkingdomarr03cuvi')
	uu <- bhl_getitembyidentifier(type='ia', value='animalkingdomarr03cuvi', as='xml')
	vv <- bhl_getitembyidentifier(type='ia', value='animalkingdomarr03cuvi', as='list')

	expect_is(tt$data, "data.frame")

	expect_is(uu, "character")
	expect_is(xmlParse(uu), "XMLInternalDocument")
	expect_is(xpathApply(xmlParse(uu), "//Year")[[1]], "XMLInternalElementNode")
	expect_is(xpathApply(xmlParse(uu), "//Year", xmlValue)[[1]][[1]], "character")

	expect_is(vv, "list")
	expect_is(vv$Status, "character")
	expect_null(vv$ErrorMessage)

  expect_equal(NCOL(tt$data), 22)
  expect_equal(length(uu), 1)
  expect_equal(length(xmlParse(uu)), 1)
  expect_equal(length(vv$Status), 1)
  expect_equal(length(vv), 3)
  expect_equal(length(vv$Result), 1)
})
