# NB2 maximum likelihood function J Hilbe and A Robinson 11Apr 2010, 10Jul 2011
ml.nb2 <- function(formula, data, offset = 0, start = NULL, verbose = FALSE) {
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  nb2X <- model.matrix(formula, data = data)
  nb2.reg.ml <- function(b.hat, X, y) {
    a.hat <- b.hat[1]
    xb.hat <- X %*% b.hat[-1] + offset
    mu.hat <- exp(xb.hat)
    r.hat <-  1 / a.hat    
    sum(dnbinom(y,
                size = r.hat,
                mu = mu.hat,
                log = TRUE))
    }
  if (is.null(start))
    start <- c(0.5, -1, rep(0, ncol(nb2X) - 1))
  fit <- optim(start,           
               nb2.reg.ml,
               X = nb2X,
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
  rownames(results) <- c("alpha", colnames(nb2X))
  results <- results[c(2:nrow(results), 1),]
  return(results)
}


