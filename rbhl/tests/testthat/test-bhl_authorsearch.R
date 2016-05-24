# tests for bhl_authorsearch fxn in rbhl
context("bhl_authorsearch")

test_that("bhl_authorsearch returns", {
	skip_on_cran()

	tt <- bhl_authorsearch(name='dimmock')
	uu <- bhl_authorsearch(name='dimmock', as = "list")

	# correct classes
	expect_is(tt$data, "data.frame")
	expect_is(tt$data$CreatorID, "integer")

	expect_is(uu, "list")
	expect_is(uu$Status, "character")
	expect_is(uu$Result, "list")
	expect_is(uu$Result[[1]]$Dates, "character")

  expect_is(bhl_authorsearch(name='dimmock', as="json"), "character")
	expect_is(bhl_authorsearch(name='dimmock', as="xml"), "character")

	# correct dimensions
  expect_equal(NCOL(tt$data), 12)
  expect_equal(length(uu$Status), 1)
})
