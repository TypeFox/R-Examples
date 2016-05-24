context("Testing W_gamma \n")

set.seed(12313)
gamma.v <- c(-0.1, 0, 0.1, 0.2)
random.data <- rnorm(n = 100)

test_that("specific identities for W_gamma", {
  
  for (gg in gamma.v) {
    expect_equal(W_gamma(0, gamma = gg), 0,
                 info = "0 at 0")
    if (gg < 0) {
      expect_true(is.nan(W_gamma(Inf, gamma = gg)),
                  info = paste("NaN at Inf and gamma = ", gg, " < 0"))
    } else {
      expect_equal(W_gamma(Inf, gamma = gg), Inf,
                   info = paste("Inf at Inf and gamma = ", gg, " >= 0"))
    }
    expect_equal(W_gamma(c(0.1, 0.2), gamma = gg),
                 -W_gamma(-c(0.1, 0.2), gamma = -gg),
                 info = "asymmetric with flipping gamma")
  }
  expect_identical(W_gamma(random.data, gamma = 0), random.data,
                   info = "identity if gamma = 0")
})

test_that("W_gamma input must be length 1", {
  
  # by default it's gamma = 0
  expect_equal(W_gamma(random.data, gamma = 0), W_gamma(random.data),
               info = "default is gamma = 0")
  expect_error(W(0, gamma = c(0, 1)),
               info = "gamma must be length one")
})

test_that("input to W_gamma must be numeric", {
  for (vv in list("foo", list(a = 1, b = 2))) {
    expect_error(W_gamma(vv, gamma = 0.2),
                 info = paste("input can't be of type ", class(vv)[1]))
  }
})

test_that("W_gamma is inverse of H_gamma", {
  for (gg in gamma.v) {
    if (gg < 0) {
      next
    }
    expect_equal(W_gamma(H_gamma(random.data, gamma = gg), gamma = gg), 
                 random.data,
                 info = paste0("W_gamma is inverse of H for gamma = ", gg))
  }
})

test_that("W_gamma treats branch correctly", {
  cat("non-principal branch is less than principal")
  expect_lt(W_gamma(-1, gamma = 0.1, branch = -1),
            W_gamma(-1, gamma = 0.1, branch = 0))
  
  cat("non-principal branch is greater than principal")
  expect_gt(W_gamma(1, gamma = -0.1, branch = -1),
            W_gamma(1, gamma = -0.1, branch = 0))
})

context("Testing derivatives")

test_that("deriv_W_gamma is correct", {
  eps <- 1e-7
  
  for (gg in gamma.v) {
    expect_equal(deriv_W_gamma(random.data, gamma = gg),
                 (W_gamma(random.data + eps, gamma = gg) - 
                    W_gamma(random.data - eps, gamma = gg)) / (2 * eps),
                 tol = 1e-2,
                 info = paste("derivative actually is the derivative ",
                              "(approximately) for gamma = ", gg))
    expect_equal(deriv_W_gamma(c(1, 10), gamma = gg),
                 deriv_W_gamma(-c(1, 10), gamma = -gg),
                 info = paste("asymmetric for gamma = ", gg))
  }
})