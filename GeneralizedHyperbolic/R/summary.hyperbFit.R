### Calculate summary for hyperbFit object
### Only applies when hessian is asked for
###
### DJS 11/08/06
summary.hyperbFit <- function(object, hessian = FALSE,
                              hessianMethod = c("exact", "tsHessian"), ...) {

  if (! "hyperbFit" %in% class(object))
    stop("Object must belong to class hyperbFit")
  obs <- object$obs
  param <- object$param


  if(hessian == TRUE) {
    hessian <- hyperbHessian(obs, param, hessianMethod = hessianMethod,
                             whichParam = 2)
    object$hessian <- hessian
    object$hessianMethod <- hessianMethod
  }

  if (!is.null(object$hessian)) {
    varcov <- solve(-object$hessian)
    par.ses <- sqrt(diag(varcov))
    object$sds <- par.ses
  }

  class(object) <- "summary.hyperbFit"
  return(object)
} ## End of summary.hyperbFit

### Print summary
print.summary.hyperbFit <-
  function(x, digits = max(3, getOption("digits") - 3), ...) {

  if (class(x) != "summary.hyperbFit")
    stop("Object must belong to class summary.hyperbFit")

  cat("\nData:     ", x$obsName, "\n")


  if (!is.null(x$hessian)) {
    cat ("Hessian: ", x$hessianMethod, "\n")
    print.default(x$hessian)
  }
  cat("Parameter estimates:\n")
  if (is.null(x$sds)) {
    print.default(format(x$param, digits = digits),
                  print.gap = 2, quote = FALSE)
  } else {
    ans <- format(rbind(x$param, x$sds), digits = digits)
    ans[1, ] <- sapply(ans[1, ], function(obs) paste("", obs))
    ans[2, ] <- sapply(ans[2, ],
                       function(obs) paste("(", obs, ")", sep = ""))
    dn <- dimnames(ans)
    dn[[1]] <- rep("", 2)
    dn[[2]] <-
        paste(substring("      ", 1,
                        (nchar(ans[2, ]) - nchar(dn[[2]])) %/% 2), dn[[2]])
    dn[[2]] <- paste(dn[[2]],
                     substring("      ", 1,
                               (nchar(ans[2, ]) - nchar(dn[[2]])) %/% 2))
    dimnames(ans) <- dn
    print.default(ans, print.gap = 2, quote = FALSE)
  }

  cat("Likelihood:        ", x$maxLik, "\n")
  cat("Method:            ", x$method, "\n")
  cat("Convergence code:  ", x$conv, "\n")
  cat("Iterations:        ", x$iter, "\n")
  invisible(x)
}
