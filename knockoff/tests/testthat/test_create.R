test_that('equicorrelated knockoffs have the right correlation structure', {
  n = 20; p = 10
  X = normc(rnorm_matrix(n,p))
  X_ko_default = knockoff.create(X, 'equicorrelated', randomize=F)
  X_ko_randomized = knockoff.create(X, 'equicorrelated', randomize=T)
  
  G = t(X) %*% X
  s = min(2*min(eigen(G)$values), 1)
  for (X_ko in list(X_ko_default, X_ko_randomized)) {
    expect_equal(t(X_ko) %*% X_ko, G)
    expect_equal(t(X) %*% X_ko, G - diag(s,p,p))
  }
})

# Test case from Weijie Su.
test_that('equicorrelated knockoffs are created in numerically sensitive case', {
  if (!requireNamespace('MASS', quietly=T))
    skip('MASS not available')
  
  n = 15; p = 5
  M = matrix(0, p, p)
  diag(M) = 1
  for (i in 1:p) {
    for (j in 1:p) {
      if ((i==j+1) || (j==i+1))
        M[i,j] <- 0.6
      if ((i==j+2) || (j==i+2))
        M[i,j] <- 0.1
    }
  }
  X = with_seed(2, MASS::mvrnorm(n, mu=rep(0,p), Sigma=M))
  k = 4
  
  Z = normc(X[,-k])
  Z_ko = knockoff.create(Z, 'equicorrelated')
  expect_false(any(is.nan(Z_ko)))
})

test_that('SDP knockoffs have the right correlation structure', {
  skip_on_cran()
  if (!has_cvxpy())
    skip('CVXPY not available')
  
  n = 20; p = 10
  X = normc(rnorm_matrix(n,p))
  X_ko_default = knockoff.create(X, 'sdp', randomize=F)
  X_ko_randomized = knockoff.create(X, 'sdp', randomize=T)
  
  offdiag <- function(A) A - diag(diag(A))
  G = t(X) %*% X
  tol = 1e-4
  for (X_ko in list(X_ko_default, X_ko_randomized)) {
    expect_equal(t(X_ko) %*% X_ko, G, tolerance=tol)
    expect_equal(offdiag(t(X) %*% X_ko), offdiag(G), tolerance=tol)
    expect_true(all(diag(t(X) %*% X_ko) < 1+tol))
  }
})