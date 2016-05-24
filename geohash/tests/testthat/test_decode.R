context("Test geohash decoding")

test_that("A simple lat/lng pair will correctly decoding", {

  result <- gh_decode("ezs42")
  expect_less_than(result$lat[1], 43)
  expect_more_than(result$lat[1], 42)
})

test_that("NAs are appropriately handled when decoding", {

  result <- gh_decode(c("ezs42", NA))
  expect_true(is.na(result$lat[2]))

})
