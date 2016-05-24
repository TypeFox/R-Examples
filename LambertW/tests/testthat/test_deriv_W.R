context("Testing derviative of Lambert W function\n")
set.seed(10)
pos.data <- rexp(100, rate = 1)

neg.data.in.range <- runif(n = 100, -exp(-1), 0)

test_that("mathematical identities for the derivative W function", {
  
  expect_equal(deriv_W(rep(0, length = 10)), rep(1, length = 10))
  expect_equal(deriv_W(-exp(-1)), Inf)
  expect_true(is.na(deriv_W(-exp(-1) - 0.0001)))

  # it actually computes the derivative
  eps <- 1e-5
  expect_equal((W(pos.data + eps) - W(pos.data - eps)) / (2 * eps),
               deriv_W(pos.data), tol = 1e-4)
  expect_equal((W(neg.data.in.range + eps, branch = -1) - 
                  W(neg.data.in.range - eps, branch = -1)) / (2 * eps),
               deriv_W(neg.data.in.range, branch = -1), tol = 1e-3)
  
  # derivative of inverse = inverse of derivative of original
  # f^(-1)'(x) = 1 / f'(x)
  #expect_equal(deriv_W(pos.data), 
  #             1 / deriv_xexp(pos.data, degree = 1),
  #             tol = 1e-2)
})



test_that("mathematical identities for the log of derivative hold", {
  
  expect_equal(log_deriv_W(rep(0, length = 10)), rep(0, length = 10))
  expect_equal(log_deriv_W(-exp(-1)), Inf)
  expect_true(is.na(log_deriv_W(-exp(-1) - 0.0001)))
  
  # it actually computes the log
  for (bb in c(0, -1)) {
    expect_equal(log(deriv_W(pos.data, branch = bb)),
                 log_deriv_W(pos.data, branch = bb), tol = 1e-4,
                 info = paste("for branch", bb))
    expect_equal(log(deriv_W(neg.data.in.range, branch = bb)), 
                 log_deriv_W(neg.data.in.range, branch = bb), 
                 tol = 1e-4,
                 info = paste("for branch", bb))  
    expect_equal(log(deriv_W(c(pos.data, neg.data.in.range), branch = bb)), 
                 log_deriv_W(c(pos.data, neg.data.in.range), branch = bb), 
                 tol = 1e-4,
                 info = paste("for branch", bb))  
  }

})


test_that("mathematical identities for the derivative of log(W) hold", {
  
  expect_equal(deriv_log_W(rep(0, length = 10)), rep(Inf, length = 10))
  expect_true(is.na(deriv_log_W(0 - 0.0001)))
  
  # it actually computes the derivative
  for (bb in c(0, -1)) {
    eps <- 1e-5
    expect_equal((log_W(pos.data + eps, branch = bb) -
                     log_W(pos.data - eps, branch = bb)) / (2 * eps),
                 deriv_log_W(pos.data, branch = bb),
                 tol = 1e-4,
                 info = paste("for branch", bb))
  }
})
