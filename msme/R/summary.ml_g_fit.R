summary.ml_g_fit <- function(object, dig = 3, ...) {
  zTable <- with(object, 
              data.frame(Estimate = beta.hat,
                         SE = se.beta.hat,
                         Z = beta.hat / se.beta.hat,
                         LCL = beta.hat - 1.96 * se.beta.hat,
                         UCL = beta.hat + 1.96 * se.beta.hat))
  rownames(zTable) <- c(colnames(object$X), "Log Sigma")
  p <- length(object$fit$par)
  n <- nrow(object$X)
  df <- c(p, n, p)
  summ <- list(call = object$call,
               coefficients = zTable,
               df = df, 
               residuals = residuals(object),
               aliased = rep(FALSE, p),
               sigma = object$sigma.hat)
  class(summ) <- c("summary.ml_g_fit", "summary.lm")
  return(summ)
}
