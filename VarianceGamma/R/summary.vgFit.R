summary.vgFit <- function (object, ...)
{
    if (!class(object) == "vgFit") {
        stop("Object must belong to class vgFit")
    }
    if (!is.null(object$hessian)) {
        sds <- sqrt(diag(solve(object$hessian)))
        sds[2] <- object$param[2] * sds[2]
        sds[3] <- object$param[4] * sds[4]
        names(sds) <- c("vgC", "sigma", "theta", "nu")
        object$sds <- sds
    }
    class(object) <- "summary.vgFit"
    return(object)
}

print.summary.vgFit <- function (x, digits = max(3, getOption("digits") - 3),
                                 ...)
{
  if (!class(x) == "summary.vgFit") {
      stop("Object must belong to class summary.vgFit")
    }
    cat("\nData:     ", x$xName, "\n")
    cat("Parameter estimates:\n")
  if (is.null(x$sds)) {
        print.default(format(x$param, digits = digits),
                      print.gap = 2, quote = FALSE)
    }
    else {
        ans <- format(rbind(x$param, x$sds), digits = digits)
        ans[1, ] <- sapply(ans[1, ], function(obs) paste("", obs))
        ans[2, ] <- sapply(ans[2, ],
                           function(obs) paste("(", obs, ")", sep = ""))
        dn <- dimnames(ans)
        dn[[1]] <- rep("", 2)
        dn[[2]] <- paste(substring("      ", 1,
                                   (nchar(ans[2, ]) - nchar(dn[[2]]))%/%2),
                         dn[[2]])
        dn[[2]] <- paste(dn[[2]],
                         substring("      ", 1,
                                   (nchar(ans[2, ]) - nchar(dn[[2]]))%/%2))
        dimnames(ans) <- dn
        print.default(ans, print.gap = 2, quote = FALSE)
    }
    cat("Likelihood:        ", x$maxLik, "\n")
    cat("Method:            ", x$method, "\n")
    cat("Convergence code:  ", x$conv, "\n")
    cat("Iterations:        ", x$iter, "\n")
    invisible(x)
}

