##' @title predPref summary
##'
##' @description  summary method for predPref objects as returned by the
##' function \code{predPref}
##'
##' @param object predPref object as returned from predPref()
##' @param ... additional arguments
##' @param sig.level significance level used in hypothesis test
##' @export
summary.predPref <- function(object, ..., sig.level=0.05) {
    out <- list()
    if ( object$LRT$p.value > sig.level ) {
        out[["estimates"]] <- object$null
    } else {
        out[["estimates"]] <- object$alt
    }

    ## indices of parameter estimates
    out[["S"]] <- ncol(out[["estimates"]][["gamma"]])
    out[["s"]] <- seq_len(out[["S"]])
    out[["T"]] <- nrow(out[["estimates"]][["gamma"]])
    out[["t"]] <- seq_len(out[["T"]])
    
    out[["p.value"]] <- object$LRT$p.value
    out[["df"]] <- object$LRT$df
    out[["Lambda"]] <- object$LRT$Lambda
    out[["loglikH0"]] <- object$loglikH0
    out[["loglikH1"]] <- object$loglikH1
    out[["hypotheses"]] <- object$hypotheses
    class(out) <- "summary.predPref"
    out
}

##' @title predPref summary print
##'
##' @description printing method for the summary function for class predPref
##'
##' @param x object of class predPref
##' @param ... additional arguments
print.summary.predPref <- function(x, ...) {
    cat("\nPredator Preferences Model:\n\n")
    cat("Hypotheses:\n\tH0 =", x$hypotheses[1], "\n\tH1 =", x$hypotheses[2], "\n")
    cat("Likelihood Ratio Test:\n\t-2*(llH0-llH1) =", x$Lambda, "on", x$df, "degrees of freedom\n")
    cat("\tp-value =", x$p.value, "\n")
    cat("Parameter Estimates:\n")

    ## indices of estimated parameters 
    indices <- unlist(lapply(x[["s"]],
                             function(y) sapply(x[["t"]],
                                                function(x) paste(x, "_", y, sep=""))))
    
    if ( !is.null(x$estimates$c) ) {
        est <- cbind(c(as.vector(x$estimates$c),
                       as.vector(x$estimates$gamma)),
                     as.vector(sqrt(diag(x$estimates$var))))
        rn <- c(paste("c_", seq_along(x$estimates$c), sep=""),
                paste("gamma_", indices, sep=""))
    } else {
        if ( converged(x$estimates$lambda, x$estimates$gamma) ) {
            est <- cbind(as.vector(x$estimates$gamma),
                     as.vector(sqrt(diag(x$estimates$var))))
            rn <- paste("gamma_", indices, sep="")
        } else {
            est <- cbind(c(as.vector(x$estimates$lambda),
                           as.vector(x$estimates$gamma)),
                         as.vector(sqrt(diag(x$estimates$var))))
            rn <- c(paste("lambda_", indices, sep=""),
                    paste("gamma_", indices, sep=""))            
        }
    }

    colnames(est) <- c("estimate", "Std. Err.")
    rownames(est) <- rn
    printCoefmat(est, ...)
    invisible(x)
}
