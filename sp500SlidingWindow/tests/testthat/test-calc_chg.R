library(sp500SlidingWindow)
context("Test calc_chg")

series <- c(16.66, 16.85, 16.93, 16.98, 17.08, 17.03)
chg    <- calc_chg(series)
foo    <- c(NA, 1.0114046, 1.0047478, 1.0029533, 1.0058893, 0.9970726)

test_that("calc_chg on a simple series", {
    expect_equal(all.equal(chg, foo, tolerance = 0.01), TRUE)
})
