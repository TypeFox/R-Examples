# Code to bound the minimum eigenvalue of the optimal matrix away from 0.

h <- function(x, a) log(x) + a / x

dhdx <- function(x, a) (1 - a / x) / x

FindMinSolution <- function(a, c, tol=1e-10, max.iter=1e3, trace=1) {
  # performs newton's method to find the smaller x s.t. h(x, a) = c.
  if (h(a, a) > c)
    stop("No solution!")
  if (h(a, a) == c)
    return(a)

  # find solution that is < a:
  x.min <- 0.1 * a
  for (i in seq(max.iter)) {
    if (trace > 1)
      cat("h(", x.min, ", a) = ", h(x.min, a), fill=T)
    x.min <- x.min - (h(x.min, a) - c) / dhdx(x.min, a)
    if (x.min <= 0)
      x.min = a * tol
    else if (x.min >= a)
      x.min = a * (1 - tol)
    if ( abs(h(x.min, a) - c) < tol & h(x.min, a) >= c) {
      if (trace > 0)
        cat("Eval-bound converged in ", i, "steps to ", x.min, fill=T)
      break
    }
  }
  x.min
}

ComputeDelta <- function(S, Lambda, Sigma.tilde=NULL, trace=1) {
  # Computes a lower bound on the minimum eigenvalue of
  # the solution of the non-convex problem.
  #
  # Args:
  #  S: Sample covariance matrix
  #  Lambda: Penalty matrix
  #  Sigma.tilde: A feasible point.  Default is S.
  if (is.null(Sigma.tilde))
    Sigma.tilde <- diag(diag(S))
  p <- nrow(S)
  f.tilde <- ComputeObjective(Sigma.tilde, S, Lambda)
  minev <- min(eigen(S)$val)
  c <- f.tilde - (p - 1) * (log(minev) + 1)
  
  FindMinSolution(a=minev, c=c, trace=trace)
}
