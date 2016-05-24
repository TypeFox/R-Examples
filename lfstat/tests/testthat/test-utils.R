context("Miscellaneous Functions ")

# generating monthly sequence
x <- seq(from=as.Date("1992-01-01"), by = "months", length.out = 12)


test_that("water_year() interprets the origin argument correctly", {

  # specify origin by an integer number
  i <- sample(1:12, size = 1)
  expect_equal(water_year(x, origin = i, assign = "end"),
               expected = as.factor(c(rep(1992, i - 1), rep(1993, 13 - i))))

  # specify origin by a POSIX object
  expect_equal(water_year(x, origin = as.POSIXct("2012-06-23"), assign = "end"),
               expected = as.factor(c(rep(1992, 5), rep(1993, 7))))

  # specify origin by a month name
  expect_equal(water_year(x, origin = "may"),
               expected = as.factor(c(rep(1991, 4), rep(1992, 8))))
  expect_equal(water_year(x, origin = "Se"),
               expected = as.factor(c(rep(1992, 8), rep(1993, 4))))

  # specify origin by the abbreviation of an institution
  expect_equal(water_year(x, origin = "din"),
               expected = as.factor(c(rep(1992, 10), rep(1993, 2))))
  expect_equal(water_year(x, origin = "usgs"),
               expected = as.factor(c(rep(1992, 9), rep(1993, 3))))
})

test_that("water_year() interprets the assign argument correctly", {

  expect_equal(water_year(x, origin = 4, assign = "majority"),
               expected = as.factor(c(rep(1991, 3), rep(1992, 9))))
  expect_equal(water_year(x, origin = 4, assign = "start"),
               expected = as.factor(c(rep(1991, 3), rep(1992, 9))))

  expect_equal(water_year(x, origin = 6, assign = "majority"),
               expected = as.factor(c(rep(1991, 5), rep(1992, 7))))
  expect_equal(water_year(x, origin = 6, assign = "start"),
               expected = as.factor(c(rep(1991, 5), rep(1992, 7))))

  expect_equal(water_year(x, origin = 7, assign = "majority"),
               expected = as.factor(c(rep(1992, 6), rep(1993, 6))))
  expect_equal(water_year(x, origin = 7, assign = "end"),
               expected = as.factor(c(rep(1992, 6), rep(1993, 6))))

})

test_that("moving averages are computed correctly", {
  # if window is of length 1, output equals input
  x <- 1:10
  expect_equal(ma(x, n = 1), x)
  expect_equal(ma(x, n = 1, sides = 2), x)

  # centered window vs. past values
  expect_equal(ma(x, n = 3, sides = 2), c(NA, 2:9, NA))
  expect_equal(ma(x, n = 3), c(NA, NA, 2:9))
})
