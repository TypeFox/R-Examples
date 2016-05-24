summary.geese <- function(object, ...) {
  mean.sum <- data.frame(estimate = object$beta,
#                         nai.se = sqrt(diag(object$vbeta.naiv)),
                         san.se = sqrt(diag(object$vbeta)),
                         ajs.se = sqrt(diag(object$vbeta.ajs)),
                         j1s.se = sqrt(diag(object$vbeta.j1s)),
                         fij.se = sqrt(diag(object$vbeta.fij)))
  mean.sum$wald <- (mean.sum$estimate / mean.sum$san.se)^2
  mean.sum$p <- 1 - pchisq(mean.sum$wald, df=1)
  rownames(mean.sum) <- object$xnames
  
  corr.sum <- data.frame(estimate = object$alpha,
#                         nai.se = sqrt(diag(object$valpha.naiv)),
                         san.se = sqrt(diag(object$valpha)),
                         ajs.se = sqrt(diag(object$valpha.ajs)),
                         j1s.se = sqrt(diag(object$valpha.j1s)),
                         fij.se = sqrt(diag(object$valpha.fij)))
  corr.sum$wald <- (corr.sum$estimate / corr.sum$san.se)^2
  corr.sum$p <- 1 - pchisq(corr.sum$wald, df=1)
  if (nrow(corr.sum) > 0) rownames(corr.sum) <- object$zcor.names

  scale.sum <- data.frame(estimate = object$gamma,
                          san.se = sqrt(diag(object$vgamma)),
                          ajs.se = sqrt(diag(object$vgamma.ajs)),
                          j1s.se = sqrt(diag(object$vgamma.j1s)),
                          fij.se = sqrt(diag(object$vgamma.fij)))
  scale.sum$wald <- (scale.sum$estimate / scale.sum$san.se)^2
  scale.sum$p <- 1 - pchisq(scale.sum$wald, df=1)
  if (!is.null(object$zsca.names)) rownames(scale.sum) <- object$zsca.names

  drop <- ifelse(c(object$control$jack, object$control$j1s, object$control$fij)== 0, TRUE, FALSE)
  if (any(drop)) {
    drop <- (3:5)[drop]
    mean.sum <- mean.sum[,-drop]
    corr.sum <- corr.sum[,-drop]
    scale.sum <- scale.sum[,-drop]
  }
  
  ans <- list(mean=mean.sum, correlation=corr.sum, scale=scale.sum,
              call=object$call, model=object$model, control=object$control,
              error=object$err, clusz=object$clusz)
  class(ans) <- "summary.geese"
  ans
}

print.geese <- function(x, digits = NULL, quote = FALSE, prefix = "", ...) {
  if(is.null(digits)) digits <- options()$digits
  else options(digits = digits)
  cat("\nCall:\n")
  dput(x$call)
  cat("\nMean Model:\n")
  cat(" Mean Link:                ", x$model$mean.link, "\n")
  cat(" Variance to Mean Relation:", x$model$variance, "\n")
  cat("\n Coefficients:\n")
  print(unclass(x$beta), digits = digits)

  if (!x$model$scale.fix) {
    cat("\nScale Model:\n")
    cat(" Scale Link:               ", x$model$sca.link, "\n")
    cat("\n Estimated Scale Parameters:\n")
    print(unclass(x$gamma), digits = digits)
  }
  else cat("\nScale is fixed.\n")
  cat("\nCorrelation Model:\n")
  cat(" Correlation Structure:    ", x$model$corstr, "\n")
  if (pmatch(x$model$corstr, c("independence", "fixed"), 0) == 0) {
    cat(" Correlation Link:         ", x$model$cor.link, "\n")
    cat("\n Estimated Correlation Parameters:\n")
    print(unclass(x$alpha), digits = digits)
  }
  ##cat("\nNumber of observations : ", x$nobs, "\n")
  ##cat("\nMaximum cluster size   : ", x$max.id, "\n")

  cat("\nReturned Error Value:  ")
  cat(x$error, "\n")
  cat("Number of clusters:  ", length(x$clusz), "  Maximum cluster size:", max(x$clusz), "\n\n")
  invisible(x)
}

print.summary.geese <- function(x, digits = NULL,
                                quote = FALSE, prefix = "", ... ) {
  if(is.null(digits)) digits <- options()$digits
  else options(digits = digits)
  cat("\nCall:\n")
  dput(x$call)
  cat("\nMean Model:\n")
  cat(" Mean Link:                ", x$model$mean.link, "\n")
  cat(" Variance to Mean Relation:", x$model$variance, "\n")
  cat("\n Coefficients:\n")
  print(x$mean, digits = digits)

  if (x$model$scale.fix == FALSE) {
    cat("\nScale Model:\n")
    cat(" Scale Link:               ", x$model$sca.link, "\n")
    cat("\n Estimated Scale Parameters:\n")
    print(x$scale, digits = digits)
  }
  else cat("\nScale is fixed.\n")
  cat("\nCorrelation Model:\n")
  cat(" Correlation Structure:    ", x$model$corstr, "\n")
  if (pmatch(x$model$corstr, c("independence", "fixed"), 0) == 0) {
    cat(" Correlation Link:         ", x$model$cor.link, "\n")
    cat("\n Estimated Correlation Parameters:\n")
    print(x$corr, digits = digits)
  }

  ##cat("\nNumber of observations : ", x$nobs, "\n")
  ##cat("\nMaximum cluster size   : ", x$max.id, "\n")

  cat("\nReturned Error Value:    ")
  cat(x$error, "\n")
  cat("Number of clusters:  ", length(x$clusz), "  Maximum cluster size:", max(x$clusz), "\n\n")
  invisible(x)
}
