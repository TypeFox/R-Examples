testthat::context("Test OLC encoding")

# Uses https://github.com/google/open-location-code/blob/master/test_data/encodingTests.csv
testthat::test_that("OLC encoding works for very simple cases", {
  testthat::expect_equal(encode_olc(20.375, 2.775, 6), "7FG49Q00+")
  testthat::expect_equal(encode_olc(20.3700625, 2.7821875, 10), "7FG49QCJ+2V")
})

testthat::test_that("OLC encoding works for more-complex cases with lengths > 10",{
  testthat::expect_equal(encode_olc(20.3701125, 2.782234375, 11), "7FG49QCJ+2VX")
  testthat::expect_equal(encode_olc(20.3701135, 2.78223535156, 13), "7FG49QCJ+2VXGJ")
  testthat::expect_equal(encode_olc(47.0000625,8.0000625, 10), "8FVC2222+22")
})

testthat::test_that("Negative values can be encoded", {
  testthat::expect_equal(encode_olc(-41.2730625,174.7859375, 10), "4VCPPQGP+Q9")
})

testthat::test_that("OLC encoding throws an error with odd-numbered requested lengths below 8",{
  testthat::expect_that(encode_olc(20.375, 2.775, 5), testthat::throws_error("The length value"))
  testthat::expect_that(encode_olc(20.375, 2.775, 7), testthat::throws_error("The length value"))
})

testthat::test_that("OLC encoding has NA support", {
  testthat::expect_that(encode_olc(20.364, NA, 5), is.na)
  testthat::expect_that(encode_olc(20.364, 2.775, NA), is.na)
  testthat::expect_that(encode_olc(NA, 2.775, 5), is.na)
})
