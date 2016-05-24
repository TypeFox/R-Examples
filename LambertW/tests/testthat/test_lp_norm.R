context("Testing Lp norm\n")

test_that("lp norm is correct for real values", {
  # Pythagoras
  expect_equal(lp_norm(c(3, 4), p = 2), 5)
  expect_equal(lp_norm(c(3, 4), p = 1), 7)
  expect_equal(lp_norm(c(3, 4), p = 0), 2)
  expect_equal(lp_norm(c(3, 4), p = Inf), 4)
  
  # 1.5 norm is between 1 and 2
  expect_lt(lp_norm(c(3, 4), p = 1.5),
            lp_norm(c(3, 4), p = 1))
  expect_gt(lp_norm(c(3, 4), p = 1.5),
            lp_norm(c(3, 4), p = 2))
  
  expect_error(lp_norm(10, -1))
})

nonneg.norms <- c(0, 1, 2, rexp(5), 10, Inf)

test_that("correct for complex values", {
  kComplexVec <- exp(1i * runif(20, -pi, pi))
  expect_equal(sapply(kComplexVec, lp_norm), rep(1, length(kComplexVec)))
  
  aa <- rnorm(10)
  for (pp in nonneg.norms) {
    info.text <- paste0("L", pp, " norm")
    expect_equal(lp_norm(aa, p = pp), lp_norm(aa + 0*1i, p = pp),
                 info = info.text)
    expect_equal(lp_norm(aa, p = pp), lp_norm(0 + 1i*aa, p = pp),
                 info = info.text)
  }
})

test_that("0 vector has 0 norm for any p", {
  for (pp in nonneg.norms) {
    expect_equal(lp_norm(rep(0, 10), p = pp), 0)
  }
})


test_that("norm of scaled vector is equal to scale times norm of vector", {
  xx <- rnorm(100)
  all.scales <- seq(-2, 2, length = 11)
  for (ss in all.scales) {
    for (pp in nonneg.norms) {
      if (pp > 0) { 
        expect_equal(lp_norm(xx * ss, p = pp),
                     lp_norm(xx, p = pp) * abs(ss),
                     info = paste("p = ", pp, "; s =", ss))
      }
    }
  }
})

