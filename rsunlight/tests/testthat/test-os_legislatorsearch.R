context("os_legislatorsearch")

test_that("returns the correct", {
  skip_on_cran()

  a <- os_legislatorsearch(state = 'ca', party = 'democratic')
  b <- os_legislatorsearch(state = 'dc', chamber = 'upper')

  # classes
  expect_is(a, "data.frame")
  expect_is(b, "data.frame")
  expect_is(a$state, "character")
  expect_is(a$`+district_offices`, "list")

  # values
  expect_equal(unique(a$state), "ca")
  expect_equal(unique(a$party), "Democratic")
  expect_equal(unique(b$state), "dc")
  expect_equal(unique(b$chamber), "upper")
})

test_that("vectorizing works", {
  skip_on_cran()

  res <- os_legislatorsearch(state = c('dc', 'or'), chamber = 'upper')
  expect_is(res, "data.frame")
  expect_equal(unique(res$state), c('dc', 'or'))
})

test_that("fails well", {
  skip_on_cran()

  expect_equal(NROW(os_legislatorsearch(state = 4)), 0)
})
