summary.simex <-
  function (object, ...)
  {
    p.names <- names(coef(object))
    est <- coef(object)
    est.table <- list()
    n <- length(resid(object))
    p <- length(p.names)
    rdf <- n - p
    if (any(names(object) == "variance.jackknife")) {
      se <- sqrt(diag(object$variance.jackknife))
      tval <- est / se
      pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
      est.table[["jackknife"]] <- cbind(est, se, tval, pval)
      dimnames(est.table[["jackknife"]]) <- list(p.names,
                                                 c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    }
    if (any(names(object) == "variance.asymptotic")) {
      se <- sqrt(diag(object$variance.asymptotic))
      tval <- est / se
      pval <- 2 * pt(abs(tval), rdf, lower.tail = FALSE)
      est.table[["asymptotic"]] <- cbind(est, se, tval, pval)
      dimnames(est.table[["asymptotic"]]) <- list(p.names,
                                                  c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    }
    ans <- list()
    class(ans) <- "summary.simex"
    ans$coefficients <- est.table
    ans$residuals <- resid(object)
    ans$call <- object$call
    ans$B <- object$B
    ans$naive.model <- object$model$call
    ans$SIMEXvariable <- object$SIMEXvariable
    ans$measurement.error <- object$measurement.error
    return(ans)
  }

