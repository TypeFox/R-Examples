context("Simconf")

test_that("Simconf", {
  data <- integration.testdata1()
  res = simconf(Q = data$Q, mu = data$mu, seed = data$seed, alpha= 0.1)
  ra = c(-7.6070830, -6.6203520, -5.6204871, -4.6204885, -3.6204885,
         -2.6204885, -1.6204885, -0.6204885,  0.3795129, 1.3796480,  2.3929170)
  rb = c(-2.3929170, -1.3796480, -0.3795129, 0.6204885, 1.6204885, 2.6204885,
         3.6204885, 4.6204885, 5.6204871, 6.6203520, 7.6070830)
  expect_equal(res$a,ra,tolerance=1e-3)
  expect_equal(res$b,rb,tolerance=1e-3)
})
