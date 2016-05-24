linprogS <- function(Sigma, e, lambda) {
  p <- nrow(Sigma)
  if (p!=ncol(Sigma)) stop("Sigma should be a square matrix!")

  f.obj <- rep(1, 2*p)
  con1 <- cbind(-Sigma, +Sigma)
  b1 <- lambda - e
  b2 <-  lambda + e
  f.con <- rbind(-diag(2*p), con1, -con1)
  f.dir <- rep("<=", 4*p)
  f.rhs <- c(rep(0,2*p), b1, b2)
  lp.out <- lp("min", f.obj, f.con, f.dir, f.rhs)
  beta <- lp.out$solution[1:p] - lp.out$solution[(p+1):(2*p)]
  if (lp.out$status == 2) warning("No feasible solution!  Try a larger lambda maybe!")
  return(beta)
}
