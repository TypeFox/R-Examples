context("Testing W_delta \n")

delta.v <- c(-0.1, 0, 0.1, 0.2)
random.data <- rnorm(n = 100)
random.data <- sort(random.data, decreasing = FALSE)
pos.data <- random.data[random.data >= 0]
neg.data <- random.data[random.data < 0]

test_that("specific identities for W_delta", {
  
  for (dd in delta.v) {
    expect_equal(W_delta(0, delta = dd), 0,
                 info = "at zero its zero for any delta")
    if (dd < 0) {
      expect_true(is.na(W_delta(Inf, delta = dd)),
                  info = "NA for negative delta and Inf input")
    } else {
      expect_equal(W_delta(Inf, delta = dd), Inf,
                   info = "W_delta(Inf, delta > 0) = Inf")
    }
    expect_equal(W_delta(pos.data, delta = dd), 
                 -W_delta(-pos.data, delta = dd),
                 info = "W_delta is symmetric around zero")

  }
  expect_identical(W_delta(random.data, delta = 0), random.data,
                   info = "identity for delta = 0")
  

})

test_that("W_delta input must be length 1", {
  
  # by default it's delta = 0
  expect_equal(W_delta(random.data, delta = 0), W_delta(random.data),
               info = "by default delta = 0")
  expect_error(W(0, delta = c(0, 1)),
               info = "delta must be length 1")
})

test_that("input to W_delta must be numeric", {
  for (vv in list("foo")) {
    expect_error(W_delta(vv, delta = 0.2))
  }
  
})

test_that("W_delta is inverse of G_delta", {
  for (dd in delta.v) {
    if (dd < 0) {
      next
    }
    expect_equal(W_delta(G_delta(random.data, delta = dd), delta = dd), 
                 random.data,
                 info = "W_delta is inverse of G_delta")
    # with alpha too
    expect_equal(W_delta_alpha(G_delta_alpha(random.data, delta = dd,
                                             alpha = 0.232),
                               delta = dd, alpha = 0.232),
                 random.data,
                 info = "W_delta_alpha is inverse of G_delta_alpha")
  }
})

test_that("double tail functions work correctly", {

  alpha.c <- 0.232
  for (dd in delta.v) {
    if (dd < 0) {
      next
    }
    w.r2 <- W_2delta(random.data, delta = c(dd, dd + 0.1))
    w.r2.split <- c(W_delta(neg.data, delta = dd), 
                    W_delta(pos.data, delta = dd + 0.1))
    expect_equal(w.r2, w.r2.split)
    
    # with alpha too
    w.r2.a <- W_2delta_2alpha(random.data, 
                              delta = c(dd, dd + 0.1),
                              alpha = c(alpha.c, alpha.c + 0.1))
    w.r2.a.split <- c(W_delta_alpha(neg.data, delta = dd, alpha = alpha.c), 
                      W_delta_alpha(pos.data, delta = dd + 0.1, 
                                    alpha = alpha.c + 0.1))
    
    expect_equal(w.r2.a, w.r2.a.split)
  }
})


test_that("deriv_W_delta is correct", {
  eps <- 1e-7
  
  for (dd in delta.v) {
    if (dd < 0) {
      next
    }
    expect_equal(deriv_W_delta(random.data, delta = dd),
                 (W_delta(random.data + eps, delta = dd) - 
                    W_delta(random.data - eps, delta = dd)) / (2 * eps),
                 tol = 1e-4,
                 info = "derivative actually computes derivative")
    expect_equal(deriv_W_delta(pos.data, delta = dd),
                 deriv_W_delta(-pos.data, delta = dd),
                 info = "derivative is symmetric")
  }
})




test_that("deriv_W_delta_alpha is correct", {
  eps <- 1e-7
  alpha.c <- 0.2312
  for (dd in delta.v) {
    if (dd < 0) {
      next
    }
    expect_equal(deriv_W_delta_alpha(random.data, delta = dd, alpha = alpha.c),
                 (W_delta_alpha(random.data + eps, delta = dd, alpha = alpha.c) - 
                    W_delta_alpha(random.data - eps, delta = dd, alpha = alpha.c)) / (2 * eps),
                 tol = 1e-4)
    expect_equal(deriv_W_delta_alpha(pos.data, delta = dd, alpha = alpha.c),
                 deriv_W_delta_alpha(-pos.data, delta = dd, alpha = alpha.c))
  }
})