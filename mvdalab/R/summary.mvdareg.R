summary.mvdareg <- function(object, ncomp = object$ncomp, digits = 3, ...) {
  x <- object
  nobj <- nrow(x$scores)
  nresp <- length(x$Ymeans)
  yvarnames <- "y" 
  if(object$val.method == "none") {
    z <- x
    ans <- z[c("call", "terms")]
    class(ans) <- "summary.mvda"
    ans$coefficients <- matrix(NA, 0L, 3L)
    dimnames(ans$coefficients) <- list(NULL, c("Estimate", "Bootstrap Error", "'t value'"))
    ans$sigma <- "N/A" 
    se <- 0 
    est <- z$coefficients[, ncomp]
    Order <- names(sort(abs(est), decreasing = T))
    est <- est[Order]
    tval <- 0 
    ans <- z[c("call", "terms")]
    ans$coefficients <- cbind(est, se, tval)
    dimnames(ans$coefficients) <- list(Order, 
                                       c("Estimate", "Bootstrap Error", "'t value'"))
    ans$sigma <- z$validation$RMSPRESS[ncomp] 
    class(ans) <- "summary.mvda"
    cat("Call:\n\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print((ans$coefficients))
    cat("\nFit Summary: \n\nNumber of objects =", nobj)
    cat("\nNumber of predictor variables =", length(attr(object$terms, "term.labels"))) 
    cat("\nMethod:", x$method)
    cat("\nDesign Matrix for Factors =", object$contrasts)
    cat("\nNumber of components considered =", ncomp)
    R2X <- sapply(1:ncomp, function(x) {
      1 - sum(diag(crossprod(as.matrix(object$Xdata) - 
                               (object$scores[ , 1:x] %*% t(object$loadings[ , 1:x]))))) / 
        sum(diag(crossprod(as.matrix(object$Xdata))))
    })
    R2Y <- 1 - sapply(1:ncomp, function(x) crossprod(object$Yactual - object$iPreds[, x])) / 
      crossprod(object$Yactual - mean(object$Yactual))
    cat("\nR2X =", round(R2X, digits))
    cat("\nR2Y =", round(R2Y, digits))
    cat("\nNo Cross-Validation")
  } else if(object$val.method == "loo") {
    z <- x
    ans <- z[c("call", "terms")]
    class(ans) <- "summary.mvda"
    ans$coefficients <- matrix(NA, 0L, 3L)
    dimnames(ans$coefficients) <- list(NULL, c("Estimate", "Bootstrap Error", "'t value'"))
    ans$sigma <- "N/A" 
    se <- 0 
    est <- z$coefficients[, ncomp]
    Order <- names(sort(abs(est), decreasing = T))
    est <- est[Order]
    tval <- 0 
    ans <- z[c("call", "terms")]
    ans$coefficients <- cbind(est, se, tval)
    dimnames(ans$coefficients) <- list(Order, 
                                       c("Estimate", "Bootstrap Error", "'t value'"))
    ans$sigma <- z$validation$RMSPRESS[ncomp] #sqrt(resvar)
    class(ans) <- "summary.mvda"
    cat("Call:\n\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print((ans$coefficients))
    cat("\nFit Summary: \n\nNumber of objects =", nobj)
    cat("\nNumber of predictor variables =", length(attr(object$terms, "term.labels"))) 
    cat("\nMethod:", x$method)
    cat("\nDesign Matrix for Factors =", object$contrasts)
    cat("\nNumber of components considered =", ncomp)
    R2X <- sapply(1:ncomp, function(x) {
      1 - sum(diag(crossprod(as.matrix(object$Xdata) - 
                               (object$scores[ , 1:x] %*% t(object$loadings[ , 1:x]))))) / 
        sum(diag(crossprod(as.matrix(object$Xdata))))
    })
    R2Y <- 1 - sapply(1:ncomp, function(x) crossprod(object$Yactual - object$iPreds[, x])) / 
      crossprod(object$Yactual - mean(object$Yactual))
    cat("\nR2X =", round(R2X, digits))
    cat("\nR2Y =", round(R2Y, digits))    
    cat("\nNo. of LOO samples = ", nobj)
    cat("\nNumber of components considered\nin above parameter estimates =", ncomp)
    cat("\nCV R2 (per component) =", round(x$validation$cvR2[1:ncomp], digits))
    cat("\nPRESS (per component) =", round(x$validation$PRESS[1:ncomp], digits))
    cat("\nPredictive Squared Error (per component) =", round(x$validation$MSPRESS[1:ncomp], digits))
    cat("\nstandard deviation error in prediction (per component) =", round(x$validation$RMSPRESS[1:ncomp], digits))
    cat("\nLOO Cross-Validation\nIf you want standard errors you need to run 'oob'")
  } else {
    z <- x
    ans <- z[c("call", "terms")]
    class(ans) <- "summary.mvda"
    ans$coefficients <- matrix(NA, 0L, 3L)
    dimnames(ans$coefficients) <- list(NULL, c("Estimate", "Bootstrap Error", "'t value'"))
    ans$sigma <- z$validation$RMSPRESS[ncomp]
    coefficient.boots <- do.call("rbind", z$validation$coefficients)
    coefficient.boots <- as.matrix(coefficient.boots[, 1:ncomp])
    cols.for.se <- z$validation$boot.means
    se <- cols.for.se[[ncomp]][, c(4, 3)]; names(se) <- c("variables", "bootsterror")
    est <- as.matrix(z$coefficients[, ncomp])
    est <- data.frame(variables = row.names(est), estimate = est); row.names(est) <- NULL
    est$variables <- as.factor(est$variables)
    me <- cols.for.se[[ncomp]][, c(4, 2)]
    se$variables <- as.factor(est$variables)
    me$variables <- as.factor(est$variables)
    combined <- merge(se, est, by = "variables")
    combined <- merge(combined, me, by = "variables")
    tval <- combined$estimate/combined$bootsterror
    bias <- combined$boot.mean - combined$estimate
    bias.se <- bias / combined$bootsterror
    ans$coefficients <- cbind(combined$estimate, 
                              combined$bootsterror, tval, 
                              bias, bias.se)
     dimnames(ans$coefficients) <- list(combined$variables, 
                                       c("Estimate", "Bootstrap Error", "'t value'", 
                                         "bias", "'bias t value'"))
    ans$sigma <- z$validation$RMSPRESS[ncomp] 
    Order <- names(sort(abs((ans$coefficients[, 1])), decreasing = T))
    ans$coefficients <- ans$coefficients[Order, ]
    class(ans) <- "summary.mvda"
    cat("Call:\n\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print((ans$coefficients))
    cat("\nFit Summary: \n\nNumber of objects =", nobj)
    cat("\nNumber of predictor variables =", length(attr(object$terms, "term.labels"))) 
    cat("\nMethod:", x$method)
    cat("\nDesign Matrix for Factors =", object$contrasts)
    cat("\nNo. of bootstrap samples = ", x$validation$bootstraps)
    cat("\nNumber of components considered\nin above parameter estimates =", ncomp)
    R2X <- sapply(1:ncomp, function(x) {
        1 - sum(diag(crossprod(as.matrix(object$Xdata) - 
                (object$scores[ , 1:x] %*% t(object$loadings[ , 1:x]))))) / 
        sum(diag(crossprod(as.matrix(object$Xdata))))
          })
    R2Y <- 1 - sapply(1:ncomp, function(x) crossprod(object$Yactual - object$iPreds[, x])) / 
          crossprod(object$Yactual - mean(object$Yactual))
    cat("\nR2X =", round(R2X, digits))
    cat("\nR2Y =", round(R2Y, digits))
    cat("\nOut-of-Bag R2 (per component) =", round(x$validation$cvR2[1:ncomp], digits))
    cat("\nOut-of-Bag PRESS (per component) =", round(x$validation$PRESS[1:ncomp], digits))
    cat("\nOut-of-Bag MSPRESS.632 (per component) =", round(x$validation$MSPRESS.632[1:ncomp], digits))
    cat("\nOut-of-Bag RMSPRESS.632 (per component) =", round(x$validation$RMSPRESS.632[1:ncomp], digits))
  }
}