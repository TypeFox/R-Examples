context("Test geohash neighbour generation")

test_that("All neighbours can be found for a simple hash", {

  result <- gh_neighbours("ezs42")
  expect_that(result$north[1], equals("ezs48"))
  expect_that(result$northeast[1], equals("ezs49"))
  expect_that(result$east[1], equals("ezs43"))
  expect_that(result$southeast[1], equals("ezs41"))
  expect_that(result$south[1], equals("ezs40"))
  expect_that(result$southwest[1], equals("ezefp"))
  expect_that(result$west[1], equals("ezefr"))
  expect_that(result$northwest[1], equals("ezefx"))
})

test_that("NAs are appropriately handled with neighbouring", {

  result <- gh_neighbours(c("ezs42", NA))
  expect_true(is.na(result$north[2]))

})

test_that("Individual neighbour extraction works", {
  hash <- "ezs42"
  expect_that(north(hash), equals("ezs48"))
  expect_that(northeast(hash), equals("ezs49"))
  expect_that(east(hash), equals("ezs43"))
  expect_that(southeast(hash), equals("ezs41"))
  expect_that(south(hash), equals("ezs40"))
  expect_that(southwest(hash), equals("ezefp"))
  expect_that(west(hash), equals("ezefr"))
  expect_that(northwest(hash), equals("ezefx"))
})
