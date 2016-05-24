# generate data
n = 100
coords = matrix(runif(n*2), nrow = n, ncol = 2)

# create covariance matrix
d = as.matrix(dist(coords))
V = 3*exp(-d/2) + 0.1*diag(n)

d1 = decomp.cov(V, "chol")
d2 = decomp.cov(V, "eigen")
d3 = decomp.cov(V, "svd")

range(V - tcrossprod(d1))
range(V - tcrossprod(d2))
range(V - tcrossprod(d3))

test_that("decomp.cov takes valid arguments",{
  # "V should be a matrix or Matrix"
  expect_that(decomp.cov(as.data.frame(V), "blah"), throws_error())
  # "V must be a square numeric matrix"
  expect_that(decomp.cov(V[-1,], "blah"), throws_error())
  # "method must be 'chol', 'eigen', or 'svd'"
  expect_that(decomp.cov(V, "blah"), throws_error())
})

test_that("decomp.cov is accurate for Matrix class matrices",{
  expect_true(max(abs(range(V - tcrossprod(d1)))) < 1e-10)
  expect_true(max(abs(range(V - tcrossprod(d2)))) < 1e-10)
  expect_true(max(abs(range(V - tcrossprod(d3)))) < 1e-10)
})

V = Matrix(V)
d1 = decomp.cov(V, "chol")
d2 = decomp.cov(V, "eigen")
d3 = decomp.cov(V, "svd")

test_that("decomp.cov is accurate for Matrix matrices",{
  expect_true(max(abs(range(V - tcrossprod(d1)))) < 1e-10)
  expect_true(max(abs(range(V - tcrossprod(d2)))) < 1e-10)
  expect_true(max(abs(range(V - tcrossprod(d3)))) < 1e-10)
})


