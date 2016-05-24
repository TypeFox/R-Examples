context("redenominate")

.dat <- new.env()

currency(c("USD", "CCY"))
exchange_rate("USDCCY")
stock("xus", "USD")
stock("xcy", "CCY")

with(.dat, {
  xus <- xts(1:5, as.Date("2010-01-01")-5:1) 
  xcy <- xts(c(40, 45), as.Date("2010-01-01")-2:1)
  USDCCY <- xts(c(10, 11, 9, 10, 12), as.Date("2010-01-01")-5:1)
})

test_that("convert TO USD", {
  expect_equal(as.numeric(redenominate("xcy", "USD", env=.dat)), c(4.00, 3.75))
})
test_that("convert FROM USD", {
  expect_equal(as.numeric(redenominate("xus", "CCY", env=.dat)), c(10, 22, 27, 40, 60))
})
test_that("xts instead of name", {
  expect_equal(as.numeric(redenominate(.dat$xcy, "USD", "CCY", env=.dat)), c(4.00, 3.75))
})
test_that("uses currency from instrument if old_base is missing", {
  expect_equal(as.numeric(redenominate("xus", "CCY", env=.dat)), c(10, 22, 27, 40, 60))
})
test_that("provided old_base overrides instrument", {
  expect_equal(as.numeric(redenominate("xus", "CCY", "USD", env=.dat)), c(10, 22, 27, 40, 60))
})
test_that("no instrument defined", {
  rm_instruments("xus")
  expect_equal(as.numeric(redenominate("xus", "CCY", "USD", env=.dat)), c(10, 22, 27, 40, 60))
  expect_error(redenominate("xus", "CCY", env=.dat), "old_base is not provided")
})
test_that("inverts FX if necessary", {
  with.usdccy <- redenominate("xus", "CCY", "USD", env=.dat)
  with(.dat, {
    CCYUSD <- 1 / USDCCY
    rm(USDCCY)
  })
  expect_equal(with.usdccy, redenominate("xus", "CCY", "USD", env=.dat))
})

# TODO: test mixing intraday and daily data