context("Testing H and H_gamma function\n")

gamma.v <- seq(-0.2, 0.3, by = .1)
random.data <- rnorm(100)

test_that("specific mathematical identities for the xexp function", {
  
  expect_equal(xexp(0), 0)
  expect_equal(xexp(1), exp(1))
  expect_equal(xexp(c(0, 1)), c(0, exp(1)))
  
  expect_equal(xexp(-1), - exp(-1))
  
  # vectorized
  zero.mat <- matrix(0, ncol = 5, nrow = 4)
  expect_identical(xexp(zero.mat), zero.mat)
})

test_that("H_gamma", {
  for (gg in gamma.v) {
    expect_equal(H_gamma(0, gamma = gg), 0)
    # identities
    if (gg != 0) {
      expect_equal(H_gamma(random.data, gamma = gg), 
                   xexp(gg * random.data) / gg)
    } else {
      expect_equal(H_gamma(random.data, gamma = gg), 
                   random.data)
    }
  }
})

test_that("derivative for xexp function", {
  
  expect_equal(deriv_xexp(0), 1)
  expect_equal(deriv_xexp(random.data), 
               exp(random.data) * random.data + exp(random.data))
  
  # zero derivative is the actual function
  expect_equal(deriv_xexp(random.data, 0),
               xexp(random.data))
  # vectorized
  zero.mat <- matrix(0, ncol = 5, nrow = 4)
  expect_identical(xexp(zero.mat), zero.mat)
  
  eps <- 1e-5
  expect_equal((xexp(random.data + eps) - xexp(random.data - eps)) / (2 * eps),
               deriv_xexp(random.data), tol = 1e-4)
  
  # second deriv
  expect_equal((deriv_xexp(random.data + eps, degree = 1) - 
                  deriv_xexp(random.data - eps, degree = 1)) / (2 * eps),
               deriv_xexp(random.data, degree = 2), tol = 1e-4)
})


