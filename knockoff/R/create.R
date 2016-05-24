#' Create knockoff variables
#' 
#' Creates knockoff variables for the original variables.
#' 
#' @param X normalized n-by-p design matrix (n >= 2p)
#' @param method either 'equicorrelated' or 'sdp'
#' @param randomize whether the knockoffs are deterministic or randomized
#' @return The n-by-p knockoff matrix
#' 
#' @details To use SDP knockoffs, you must have a Python installation with 
#' CVXPY. For more information, see the vignette on SDP knockoffs:
#'
#' \code{vignette('sdp', package='knockoff')}
#' 
#' @export
knockoff.create <- function(X, method=c('equicorrelated','sdp'), randomize=F) {
  fn = switch(match.arg(method), 
              equicorrelated = create_equicorrelated,
              sdp = create_sdp)
  fn(X, randomize)
}

# Create equicorrelated knockoffs.
create_equicorrelated <- function(X, randomize) {
  # Compute SVD and U_perp.
  X.svd = decompose(X, randomize)
  
  # Set s = min(2 * smallest eigenvalue of X'X, 1), so that all the correlations
  # have the same value 1-s.
  if (any(X.svd$d <= 1e-5 * max(X.svd$d)))
    stop(paste('Data matrix is rank deficient.',
               'Equicorrelated knockoffs will have no power.'))
  lambda_min = min(X.svd$d)^2
  s = min(2*lambda_min, 1)
  
  # Construct the knockoff according to Equation 1.4.
  s_diff = pmax(0, 2*s - (s/X.svd$d)^2) # can be negative due to numerical error
  X_ko = (X.svd$u %*diag% (X.svd$d - s / X.svd$d) +
          X.svd$u_perp %*diag% sqrt(s_diff)) %*% t(X.svd$v)
}

# Create SDP knockoffs.
create_sdp <- function(X, randomize) {
  if (!has_cvxpy())
    stop('To use SDP knockoffs, you must have Python and the CVXPY library.')
  
  # Compute SVD and U_perp.
  X.svd = decompose(X, randomize)
  
  # Check for rank deficiency.
  tol = 1e-5
  d = X.svd$d
  d_inv = 1 / d
  d_zeros = d <= tol*max(d)
  if (any(d_zeros)) {
    warning(paste('Data matrix is rank deficient.',
                  'Model is not identifiable, but proceeding with SDP knockoffs'))
    d_inv[d_zeros] = 0
  }
  
  # Compute the Gram matrix and its (pseudo)inverse.
  G = (X.svd$v %*diag% d^2) %*% t(X.svd$v)
  G_inv = (X.svd$v %*diag% d_inv^2) %*% t(X.svd$v)
  
  # Optimize the parameter s of Equation 1.3 using SDP.
  s = solve_sdp(G)
  s[s <= tol] = 0
  
  # Construct the knockoff according to Equation 1.4:
  C.svd = canonical_svd(2*diag(s) - (s %diag*% G_inv %*diag% s))
  X_ko = X - (X %*% G_inv %*diag% s) + 
    (X.svd$u_perp %*diag% sqrt(pmax(0, C.svd$d))) %*% t(C.svd$v)
}

# Compute the SVD of X and construct an orthogonal matrix U_perp such that
# U_perp * U = 0.
decompose <- function(X, randomize) {
  n = nrow(X); p = ncol(X)
  stopifnot(n >= 2*p)
  
  result = canonical_svd(X)
  Q = qr.Q(qr(cbind(result$u, matrix(0,n,p))))
  u_perp = Q[,(p+1):(2*p)]
  if (randomize) {
      Q = qr.Q(qr(rnorm_matrix(p,p)))
      u_perp = u_perp %*% Q
  }
  result$u_perp = u_perp
  result
}

# Solves the semidefinite programming problem:
#
#   maximize    sum(s)
#   subject to  0 <= s <= 1
#               2G - diag(s) >= 0
#
# Because R lacks a decent SDP library, we call out to Python. It would be nice
# to use an R-to-Python FFI, but there is only one maintained package (rPython)
# and it suffers from several defects:
#
#   1) No MS Windows support
#
#   2) The Python interpreter is fixed at compile time. Consequently any OS X
#      user not using the system Python (i.e., every OS X user) will have to
#      compile rPython from scratch (which is not as easy as one might expect).
#
# So we call the Python interpreter directly, using JSON as the data 
# serialization format.
solve_sdp <- function(G) {
  source_file = system.file('python', 'solve_sdp.py', package='knockoff')
  out_file = tempfile(pattern='knockoff', fileext='json')
  on.exit({
    if (file.exists(out_file))
      file.remove(out_file)
  })
  
  G.json = toJSON(G, collapse=' ')
  status = system2('python',
                   paste(shQuote(source_file), '-', shQuote(out_file)),
                   input = G.json)
  if (status != 0)
    stop('Error calling Python SDP solver.')
  s.json = paste(readLines(out_file, warn=F))
  s = fromJSON(s.json)
  return(s)
}