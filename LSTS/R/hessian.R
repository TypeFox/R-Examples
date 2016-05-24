hessian = function(f, x0, ...){
  n = length(x0)
  grad = rep(0, n)
  mdelta = matrix(0., nrow = n, ncol = n)
  hess = mdelta
  delta = 0.0001 * (sign(x0) + (x0 == 0.)) * pmax(abs(x0), 0.01)
  diag(mdelta) = delta
  f0 = f(x0, ...)
  for(i in 1.:n) {
    grad[i] <- f(x0 + mdelta[, i], ...)
  }
  for(i in 1.:n) {
    for(j in i:n) {
      hess[i, j] <- f(x0 + mdelta[, i] + mdelta[, j], ...) - grad[i] - grad[j]
      hess[j, i] <- hess[i, j]
    }
  }
  (hess + f0)/outer(delta, delta, "*")
}