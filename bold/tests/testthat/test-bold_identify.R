context("bold_identify")

seq <- sequences$seq1

test_that("bold_identify works as expected", {
  skip_on_cran()

  aa <- bold_identify(seq)
  expect_is(aa, 'list')
  expect_is(aa[[1]], 'data.frame')
  expect_is(aa[[1]]$ID, 'character')
})

test_that("bold_identify db param works as expected", {
  skip_on_cran()

  aa <- bold_identify(seq, db = 'COX1_SPECIES')
  expect_is(aa, 'list')
  expect_is(aa[[1]], 'data.frame')
  expect_is(aa[[1]]$ID, 'character')
})

test_that("bold_identify response param works as expected", {
  skip_on_cran()

  aa <- bold_identify(seq, response = TRUE)
  expect_is(aa, "list")
  expect_is(aa[[1]], "response")
  expect_equal(aa[[1]]$status_code, 200)
})

test_that("bold_identify fails well", {
  skip_on_cran()

  expect_error(bold_identify(), "argument \"sequences\" is missing, with no default")
})
