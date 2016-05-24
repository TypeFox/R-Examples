summary.coxinterval <- function(object, conf.int = 0.95, scale = 1, ...)
{
  s <- sapply(c("call", "loglik", "iter", "censor.rate", "censor"),
              function(x) object[[x]], simplify = FALSE)
  f <- function(fit) {
    out <- sapply(c("n", "p", "m", "na.action"),
                  function(y) fit[[y]], simplify = FALSE)
    out$formula <- fit$call$formula
    if (is.null(fit$coef)) out$fit <- NULL
    else {
      est <- fit$coef
      se <- sqrt(diag(fit$var))
      out$coef.mat <- cbind(est, exp(est), se, est/se, 1 - pchisq((est/se)^2, 1))
      dimnames(out$coef) <-
        list(names(est), c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)"))
      if (conf.int) {
        lower <- (1 - conf.int)/2
        upper <- (1 + conf.int)/2
        est <- scale * est
        se <- scale * se
        out$conf <- cbind(exp(-est), exp(est), exp(est + qnorm(lower) * se),
                          exp(est + qnorm(upper) * se))
        dimnames(out$conf) <-
          list(names(est),
               c("exp(-coef)", "exp(coef)",
                 paste(round(100 * c(lower, upper), 2), "%", sep = "")))
      }
    }
    out
  }
  s <- c(s, f(object))
  temp <- c("coxph", "timereg")
  temp <- temp[which(temp %in% names(object))]
  s$rcfit <- if (length(temp) && !is.null(temp <- object[[temp]])) f(temp)
             else NULL
  class(s) <- "summary.coxinterval"
  s
}
