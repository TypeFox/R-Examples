context("Testing H_gamma \n")

set.seed(12313)
gamma.v <- c(-0.1, 0, 0.1, 0.2)
random.data <- rnorm(n = 100)

test_that("specific identities for H_gamma", {
  
  for (gg in gamma.v) {
    expect_equal(H_gamma(0, gamma = gg), 0,
                 info = "0 at 0")
    if (gg < 0) {
      expect_true(is.nan(H_gamma(Inf, gamma = gg)),
                  info = paste("NaN at Inf and gamma = ", gg, " < 0"))
    } else {
      expect_equal(H_gamma(Inf, gamma = gg), Inf,
                   info = paste("Inf at Inf and gamma = ", gg, " >= 0"))
    }
    expect_equal(H_gamma(c(0.1, 0.2), gamma = gg),
                 -H_gamma(-c(0.1, 0.2), gamma = -gg),
                 info = "asymmetric with flipping gamma")
  }
  expect_identical(H_gamma(random.data, gamma = 0), random.data,
                   info = "identity if gamma = 0")
})

test_that("H_gamma input must be length 1", {
  
  # by default it's gamma = 0
  expect_equal(H_gamma(random.data, gamma = 0), H_gamma(random.data),
               info = "default is gamma = 0")
  # xexp() = H()
  expect_error(xexp(0, gamma = c(0, 1)),
               info = "gamma must be length one")
})

test_that("input to H_gamma must be numeric", {
  for (vv in list("foo", list(a = 1, b = 2))) {
    expect_error(H_gamma(vv, gamma = 0.2),
                 info = paste("input can't be of type ", class(vv)[1]))
  }
})
