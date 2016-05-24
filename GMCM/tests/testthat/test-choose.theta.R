context("Check choose.theta")

set.seed(1)
u <- SimulateGMCMData(n = 100, m = 3, d = 2)$u

test_that("choose.theta returns proper formatted theta", {
  for (m in 1:10) {
    expect_that(is.theta(choose.theta(u, m = m)), is_true())
  }
})

test_that("choose.theta fails with no clusters", {
   expect_that(choose.theta(u, m = 0),
               throws_error("number of cluster centres must lie between"))
})
