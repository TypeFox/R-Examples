context("Integration")
test_that("Integration L", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q.chol = data$L, a = data$a, b = data$b,
                    seed = data$seed, max.threads = 1)
  expect_equal(prob1$P[1], 0.9680023, tolerance=1e-7)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance=1e-6)
})

test_that("Integration Q", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q = data$Q, a = data$a, b = data$b,
                    seed = data$seed, max.threads = 1)
  expect_equal(prob1$P[1], 0.9680023, tolerance=1e-7)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance=1e-6)
})

test_that("Integration mu", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q = data$Q, mu = data$mu, a = data$a + data$mu,
                    b = data$b + data$mu, seed = data$seed,
                    max.threads = 1)
  expect_equal(prob1$P[1], 0.9680023, tolerance=1e-7)
  expect_equal(prob1$E[1], 5.914764e-06, tolerance=1e-6)
})

test_that("Integration limit", {
  data <- integration.testdata1()
  prob1 <- gaussint(Q = data$Q, a = data$a, b = data$b,
                    seed = data$seed, lim = 0.97,
                    max.threads = 1)

  prob2 <- gaussint(Q = data$Q, a = data$a, b = data$b,
                   seed = data$seed, lim = 0.9,
                   max.threads = 1)

  expect_equal(prob1$P[1], 0.0, tolerance=1e-7)
  expect_equal(prob2$P[1],  0.9680023, tolerance=1e-6)
})


#test_that("Integration seed", {
#
# seed = 1:6
# n=10
# x <- excursions:::excursions.rand(n,seed,n.threads=1)

# r <- c(0.001009498, 0.5950038, 0.3578345, 0.2223408, 0.4668276,
#        0.3789078, 0.006938634, 0.9940359, 0.7597992, 0.8110961)

# expect_equal(x, r, tolerance=1e-4)
#})
