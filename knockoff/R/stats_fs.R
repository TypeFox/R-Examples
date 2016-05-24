#' Knockoff filter forward selection statistics
#' 
#' Computes the statistic
#'   \deqn{W_j = \max(Z_j, Z_{j+p}) \cdot \mathrm{sgn}(Z_j - Z_{j+p}),}
#' where \eqn{Z_1,\dots,Z_{2p}} give the reverse order in which the 2p
#' variables (the originals and the knockoffs) enter the forward selection 
#' model. See the Details for information about forward selection.
#' 
#' In \emph{forward selection}, the variables are chosen iteratively to maximize
#' the inner product with the residual from the previous step. The initial
#' residual is always \code{y}. In standard forward selection
#' (\code{knockoff.stats.fs}), the next residual is the remainder after
#' regressing on the selected variable; when \emph{orthogonal matching pursuit}
#' is used (\code{knockoff.stats.fs_omp}), the next residual is the remainder
#' after regressing on \emph{all} the previously selected variables.
#' 
#' @param X original design matrix
#' @param X_ko knockoff matrix
#' @param y response vector
#' @return The statistic W
#' 
#' @rdname knockoff.stat.fs
#' @export
knockoff.stat.fs <- function(X, X_ko, y) {
  stat.fs(X, X_ko, y, omp=FALSE)
}
#' @rdname knockoff.stat.fs
#' @export
knockoff.stat.fs_omp <- function(X, X_ko, y) {
  stat.fs(X, X_ko, y, omp=TRUE)
}
# Internal implementation.
stat.fs <- function(X, X_ko, y, omp) {
  p = ncol(X)
  path = fs(cbind(X, X_ko), y, omp)
  Z = 2*p + 1 - order(path)
  orig = 1:p
  pmax(Z[orig], Z[orig+p]) * sign(Z[orig] - Z[orig+p])
}

#' Forward selection
#' 
#' Perform forward variable selection with or without OMP
#' 
#' @param X matrix of predictors
#' @param y response vector
#' @param omp whether to use orthogonal matching pursuit (OMP)
#' @return vector with jth component the variable added at step j
#' 
#' @keywords internal
fs <- function(X, y, omp = FALSE) {
  n = nrow(X); p = ncol(X)
  stopifnot(n == length(y))
  
  path = rep.int(0, p)
  in_model = rep(FALSE, p)
  residual = y
  if (omp) Q = matrix(0, n, p)
  
  for (step in 1:p) {
    # Find the best variable to add among the remaining variables.
    available_vars = which(!in_model)
    products = apply(X[,!in_model,drop=F], 2,
                     function(x) abs(sum(x * residual)))
    best_var = available_vars[which.max(products)]
    path[step] = best_var
    in_model[best_var] = TRUE
    
    # Update the residual.
    x = X[,best_var]
    if (step == p) break
    if (omp) {
      for (j in seq(1, length.out=step-1))
        x = x - sum(Q[,j]*x) * Q[,j]
      q = x / sqrt(sum(x^2))
      Q[,step] = q
      residual = residual - sum(q*y) * q
    } else
      residual = residual - sum(x*residual) * x
  }
  return(path)
}