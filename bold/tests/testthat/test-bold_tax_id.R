context("bold_tax_id")

test_that("bold_tax_id returns the correct classes", {
  skip_on_cran()

  aa <- bold_tax_id(88899)
  bb <- bold_tax_id(125295)

  expect_is(aa, "data.frame")
  expect_is(bb, "data.frame")

  expect_is(aa$input, "numeric")
  expect_is(aa$taxid, "integer")
  expect_is(aa$tax_rank, "character")
})

test_that("bold_tax_id works with multiple ids passed in", {
  skip_on_cran()

  aa <- bold_tax_id(c(88899,125295))

  expect_is(aa, "data.frame")
  expect_equal(NROW(aa), 2)
})

test_that("bold_tax_id dataTypes param works as expected", {
  skip_on_cran()

  aa <- bold_tax_id(88899, dataTypes = "basic")
  bb <- bold_tax_id(88899, dataTypes = "stats")
  dd <- bold_tax_id(88899, dataTypes = "geo")
  ee <- bold_tax_id(88899, dataTypes = "sequencinglabs")

  expect_is(aa, "data.frame")
  expect_is(bb, "data.frame")
  expect_is(dd, "data.frame")
  expect_is(ee, "data.frame")

  expect_equal(NROW(aa), 1)
  expect_equal(NROW(bb), 1)
  expect_equal(NROW(dd), 1)
  expect_equal(NROW(ee), 1)

  expect_named(dd, c('input','Brazil','Mexico','Panama','Guatemala','Peru','Bolivia','Ecuador'))

  expect_gt(NCOL(bb), NCOL(aa))
  expect_gt(NCOL(ee), NCOL(aa))
  expect_gt(NCOL(bb), NCOL(ee))
})

test_that("includeTree param works as expected", {
  skip_on_cran()

  aa <- bold_tax_id(id=88899, includeTree=FALSE)
  bb <- bold_tax_id(id=88899, includeTree=TRUE)

  expect_is(aa, "data.frame")
  expect_is(bb, "data.frame")
  expect_gt(NROW(bb), NROW(aa))
})

test_that("bold_tax_id fails well", {
  skip_on_cran()

  expect_error(bold_tax_id(), "argument \"id\" is missing, with no default")
})
