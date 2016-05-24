context("splitTime")

test_that("splitTime", {
  expect_equal(splitTime(0, "years"), c(years=0, days=0, hours=0, minutes=0, seconds=0))
  expect_equal(splitTime(0, "days"), c(years=NA, days=0, hours=0, minutes=0, seconds=0))
  expect_equal(splitTime(0, "hours"), c(years=NA, days=NA, hours=0, minutes=0, seconds=0))
  expect_equal(splitTime(0, "minutes"), c(years=NA, days=NA, hours=NA, minutes=0, seconds=0))
  expect_equal(splitTime(0, "seconds"), c(years=NA, days=NA, hours=NA, minutes=NA, seconds=0))

  seconds = 2 * 365 * 24 * 60 * 60
  expect_equal(splitTime(seconds, "years"), c(years=2, days=0, hours=0, minutes=0, seconds=0))
  expect_equal(splitTime(seconds, "days"), c(years=NA, days=2 * 365, hours=0, minutes=0, seconds=0))
  expect_equal(splitTime(seconds, "hours"),
    c(years=NA, days=NA, hours=2 * 365 * 24, minutes=0, seconds=0))
  expect_equal(splitTime(seconds, "minutes"),
    c(years=NA, days=NA, hours=NA, minutes=2 * 365 * 24 * 60, seconds=0))
  expect_equal(splitTime(seconds, "seconds"),
    c(years=NA, days=NA, hours=NA, minutes=NA, seconds=seconds))

  expect_true(is.integer(splitTime(100000, "minutes")))
})
