context("Testing Lambert W function\n")

test_that("specific mathematical identities for the W function", {
  zz <- rnorm(1e3)
  
  for (bb in c(0, -1)) {
    # it actually computes the log
    expect_equal(log_W(zz, branch = bb), log(W(zz, branch = bb)))
  
    # log(W(z)) = log(z) - W(z)
    expect_equal(log_W(zz, branch = bb), log(zz) - W(zz, branch = bb))
  
    expect_equal(log_W(zz, W.z = W(zz, branch = bb)), 
                 log_W(zz, branch = bb))
  }
  # at Inf it returns Inf
  expect_equal(log_W(c(0, Inf)), c(-Inf, Inf))
  
  # at 0 it is -Inf
  expect_equal(log_W(0), -Inf)
  
  # it is NaN for < 0
  expect_true(all(is.nan(log_W(-rexp(10, 1)))))
  
})
