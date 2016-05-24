context("Testing Lambert W function\n")

test_that("specific mathematical identities for the W function", {
  
  expect_equal(W(0, branch = 0), 0)
  expect_equal(W(exp(1), branch = 0), 1)
  expect_equal(W(-exp(-1), branch = 0), -1)
  
  # coincide at z = -1
  expect_equal(W(-1, branch = -1), W(-1, branch = 0))
  for (bb in c(0, -1)) {
    # for values z < -1, both branches are NA
    expect_true(all(is.na(W(-c(2:10), branch = bb))))
  }

  expect_equal(W(Inf, branch = 0), Inf)
  
  # At the critical point z = -1, both branches coincide.
  expect_equal(W(-1, branch = -1), W(-1, branch = 0))
  
  #  see also https://en.wikipedia.org/wiki/Lambert_W_function
  pos.vals <- rexp(100)
  # log(W(x)) = log(x) - W(x)
  expect_equal(log(W(pos.vals)), log(pos.vals) - W(pos.vals))
  # W(x) * exp(W(x)) = x
  expect_equal(W(pos.vals) * exp(W(pos.vals)), pos.vals)
})

test_that("W has branch hast length 1, and is either 0 or -1", {
  
  # by default it's branch 0
  expect_equal(W(c(-100, 100)), W(c(-100, 100), branch = 0))
  expect_error(W(0, branch = -2))
  expect_error(W(0, branch = c(0, 0)))
})

test_that("input to W must be numeric", {
  for (vv in list("foo")) {
    expect_error(W(vv))
  }
  
})

test_that("W is inverse of xexp", {
  pos.vals <- rexp(100, rate = 10)
  expect_equal(W(xexp(pos.vals), branch = 0), pos.vals)
  
  neg.vals <- seq(-2, -10)
  expect_equal(W(xexp(neg.vals), branch = -1), neg.vals)
})


test_that("W(z) is asymptotically like log(z)", {
  # after ~exp(4) the ration is converging monotonically from below to 1
  z.vals <- exp(c(4:20))

  W.log.ratio <- W(z.vals) / log(z.vals)
  
  expect_true(all(W.log.ratio > 0))
  expect_true(all(W.log.ratio < 1))
  # monotonically increasing
  expect_true(all(diff(W.log.ratio) > 0))
})

test_that("W throws warning if input is NA or NaN and returns NA; Inf for Inf", {

  expect_equal(W(c(0, Inf)), c(0, Inf))
  
  expect_warning(W(c(1, NaN)))
  expect_equal(W(c(0, NaN)), c(0, NA))
  expect_equal(W(c(NA, 0, NA)), c(NA, 0, NA))

})

test_that("W is vectorized and maintins input dimension", {
  
  data.list <- list("single" = 5,
                    "vector" = cbind(1:5),
                    "matrix" = matrix(rnorm(100), ncol = 2))
  
  for (dd in data.list) {
    expect_equal(dim(dd), dim(W(dd)))
  }
})