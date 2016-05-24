library(R.matlab)

# Note: This command requires that the 'matlab' binary be on the PATH.
Matlab$startServer()

matlab <- Matlab()
if (!isOpen(matlab))
  open(matlab)

call_function <- function(matlab, fn, ...) {
  args = list(...)
  names(args) = paste0('R_arg', 1:length(args))
  do.call(setVariable, c(list(matlab), args))
  
  result_name = 'R_result'
  arg_str = paste0(names(args), collapse=', ')
  call_str = paste0(result_name, '=', fn, '(', arg_str, ');')
  evaluate(matlab, call_str)
  getVariable(matlab, result_name)[[1]]
}

test_that('SVDs in R and MATLAB are identical', {
  X = normc(rnorm_matrix(15, 10))
  X.svd = canonical_svd(X)
  
  setVariable(matlab, X=X)
  evaluate(matlab, '[U,S,V] = knockoff.private.canonicalSVD(X);')
  X.svd.matlab = getVariable(matlab, c('U','S','V'))
  
  expect_equal(X.svd$u, X.svd.matlab$U)
  expect_equal(X.svd$d, diag(X.svd.matlab$S))
  expect_equal(X.svd$v, X.svd.matlab$V)
})

test_that('equicorrelated knockoffs in R and MATLAB are identical', {
  X = normc(rnorm_matrix(20, 10))
  X_ko = knockoff.create(X, 'equicorrelated')
  X_ko.matlab = call_function(matlab, 'knockoff.create', X, 'equi')
  expect_equal(X_ko, X_ko.matlab)
})

test_that('SDP knockoffs in R and MATLAB are identical', {
  if (!has_cvxpy())
    skip('CVXPY not available')
  X = normc(rnorm_matrix(20, 10))
  X_ko = knockoff.create(X, 'sdp')
  X_ko.matlab = call_function(matlab, 'knockoff.create', X, 'sdp')
  # FIXME: MATLAB and Python are using different SDP solvers, but is this
  # precision too low?
  expect_equal(X_ko, X_ko.matlab, tol=0.01)
})

test_that('forward selection in R and MATLAB are identical', {
  n = 200; p = 100
  X = normc(rnorm_matrix(n, p))
  y = X %*% rnorm(p) + 0.1 * rnorm(n)
  
  path = fs(X, y, omp=FALSE)
  path.matlab = as.vector(call_function(
    matlab, 'knockoff.private.forwardSelection', X, y))
  expect_equal(path, path.matlab)
  
  path = fs(X, y, omp=TRUE)
  path.matlab = as.vector(call_function(
    matlab, 'knockoff.private.forwardSelectionOMP', X, y))
  expect_equal(path, path.matlab)
})

test_that('test statistics in R and MATLAB are identical', {
  prob = random_problem(30, 10)
  X = prob$X; y = prob$y
  X_ko = knockoff.create(X)
  
  expect_stats_equal <- function(stat.r, stat.matlab) {
    W = stat.r(X, X_ko, y)
    W.matlab = call_function(matlab, stat.matlab, X, X_ko, y)
    expect_equal(W, c(W.matlab))
  }
  expect_stats_equal(knockoff.stat.fs, 'knockoff.stats.forwardSelection')
  expect_stats_equal(knockoff.stat.lasso_difference,
                     'knockoff.stats.lassoDifference')
  expect_stats_equal(knockoff.stat.lasso_signed_max,
                     'knockoff.stats.lassoSignedMax')
})

close(matlab)