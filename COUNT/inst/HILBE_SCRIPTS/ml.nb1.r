# NB1 maximum likelihood function  J Hilbe 11Apr 2010
ml.nb1 <- function(formula, data, start = NULL, verbose = FALSE) {
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  nb1X <- model.matrix(formula, data = data)
  nb1.reg.ml <- function(b.hat, X, y) {
    a.hat <- b.hat[1]
    xb.hat <- X %*% b.hat[-1]
    mu.hat <- exp(xb.hat)
    r.hat <- (1/a.hat) * mu.hat    
    sum(dnbinom(y,
                size = r.hat,
                mu = mu.hat,
                log = TRUE))
    }
  if (is.null(start))
    start <- c(0.5, -1, rep(0, ncol(nb1X) - 1))
  fit <- optim(start,           
               nb1.reg.ml,
               X = nb1X,
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
  rownames(results) <- c("alpha", colnames(nb1X))
  results <- results[c(2:nrow(results), 1),]
  return(results)
}

load("c://source/medpar.RData")
medpar$type <- factor(medpar$type)
med.nb1 <- ml.nb1(los ~ hmo + white + type, data = medpar)
med.nb1





