# Poisson maximum likelihood function J Hilbe and A Robinson 11Apr 2010, 10Jul 2011
ml.pois <- function(formula, data, offset = 0, start = NULL, verbose = FALSE) {
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  poisX <- model.matrix(formula, data = data)
  pois.reg.ml <- function(b.hat, X, y) {
    xb.hat <- X %*% b.hat + offset
    mu.hat <- exp(xb.hat)
    sum(dpois(y,
              lambda = mu.hat,
              log = TRUE))
    }
  if (is.null(start))
    start <- c(-1, rep(0, ncol(poisX) - 1))
  fit <- optim(start,           
               pois.reg.ml,
               X = poisX,
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
  rownames(results) <- colnames(poisX)
  results <- results[c(2:nrow(results), 1),]
  return(results)
}
