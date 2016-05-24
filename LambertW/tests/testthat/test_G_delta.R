context("Testing G_delta, G_delta_alpha function\n")

delta.v <- seq(-0.2, 0.3, by = .1)
alpha.c <- 0.234
random.data <- rnorm(100)

# vectorized
zero.mat <- matrix(0, ncol = 5, nrow = 4)

test_that("specific mathematical identities for  G_delta", {
  
  for (dd in delta.v) {
    expect_equal(G_delta(0, delta = dd), 0)
    expect_equal(G_delta(1, delta = dd), 1 * exp(dd *1/2 * 1))

    expect_identical(dim(G_delta(zero.mat, delta = dd)), 
                     dim(zero.mat))
  }
  
  expect_equal(G_delta(random.data, delta = 0), random.data,
               info = "Identity for delta = 0")
})



test_that("specific mathematical identities for  G_delta_alpha", {
  
  for (dd in delta.v) {
    expect_equal(G_delta_alpha(0, delta = dd, alpha = alpha.c), 0)
    expect_equal(G_delta_alpha(2, delta = dd, alpha.c), 
                 2 * exp(dd * 1/2 * 2^(2 * alpha.c)))
    
    expect_identical(dim(G_delta_alpha(zero.mat, delta = dd, alpha = alpha.c)), 
                     dim(zero.mat))
  }
  
  expect_equal(G_delta_alpha(random.data, delta = 0, alpha = alpha.c), 
               random.data,
               info = "Identity for delta = 0  and any alpha")
  expect_equal(G_delta_alpha(random.data, alpha = alpha.c), 
               random.data,
               info = "Identity for delta = 0 (default value) and any alpha")
})
