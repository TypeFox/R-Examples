summary.nlFit <- function(object, ...) {

  if (! "nlFit" %in% class(object))
    stop("Object must belong to class nlFit")

  if (!is.null(object$hessian)) {
    varcov <- solve(object$hessian)
    par.ses <- sqrt(diag(varcov))
    object$sds <- par.ses
  }

  class(object) <- "summary.nlFit"
  return(object)
}


print.summary.nlFit <- function(x,
                                digits = max(3, getOption("digits") - 3),
                                ...) {

  if (class(x) != "summary.nlFit")
    stop("Object must belong to class summary.nlFit")

  cat("\nData:     ", x$obsName, "\n")
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
