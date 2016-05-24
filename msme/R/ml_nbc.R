ml.nbc <- function(formula, data, start = NULL, verbose = FALSE) {
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  nbcX <- model.matrix(formula, data = data)
  nbc.reg.ml <- function(b.hat, X, y, offset) {
    a.hat <- b.hat[1]
    xb.hat <- X %*% b.hat[-1]
    mu.hat <- 1 / ((exp(-xb.hat)-1)*a.hat)
    p.hat <- 1 / (1 + a.hat*mu.hat)
    r.hat <- 1 / a.hat
    sum(dnbinom(y,
                size = r.hat,
                prob = p.hat,
                log = TRUE))
    }
  if (is.null(start))
    start <- c(0.5, -1, rep(0, ncol(nbcX) - 1))
  fit <- optim(start,           
               nbc.reg.ml,
               X = nbcX,
               y = y,
               control = list(
                 fnscale = -1,
                 maxit = 10000),
               hessian = TRUE
               )
  if (verbose | fit$convergence > 0) print(fit)
  beta.hat <- fit$par
  se.beta.hat <- sqrt(diag(solve(-fit$hessian)))
  results <- data.frame(Estimate = beta.hat,
                        SE = se.beta.hat,
                        Z = beta.hat / se.beta.hat,
                        LCL = beta.hat - 1.96 * se.beta.hat,
                        UCL = beta.hat + 1.96 * se.beta.hat)
  rownames(results) <- c("alpha", colnames(nbcX))
  results <- results[c(2:nrow(results), 1),]
  return(results)
}

