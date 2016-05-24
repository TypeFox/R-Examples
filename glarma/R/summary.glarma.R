## Summary Method
summary.glarma <- function(object, tests = TRUE, ...) {
    sum <- list()
    sum$call <- object$call
    sum$methods <- data.frame(type = object$type, iter.method = object$method,
                              Resid.type = object$residType, row.names = "")
    sum$null.deviance <- object$null.deviance
    sum$df.null <- object$df.null
    sum$deviance <- sum(object$residuals^2)
    sum$aic <- object$aic
    sum$df.residual <- NROW(object$y) - NROW(object$delta)
    sum$phi.lags <- object$phiLags
    sum$theta.lags <- object$thetaLags
    sum$pq <- object$pq
    sum$iter <- object$iter
    sum$deviance.resid <- object$residuals
    if (tests) sum$likTests <- unclass(likTests(object))
    sum$tests <- tests
    ## sum$five.num.summary <- fivenum(object$residuals)
    ## names(sum$five.num.summary) <- c("Min", "1Q", "Median", "3Q", "Max")
    sum$coefficients1 <- glarmaModelEstimates(object)[1:ncol(object$X), ]
    if (object$type == "NegBin"){
        sum$coefficients2 <-
            glarmaModelEstimates(object)[(ncol(object$X) + 1):
                                         (ncol(object$X) +
                                          length(object$thetaLags) +
                                          length(object$phiLags)), ]
        sum$coefficients3 <-
            glarmaModelEstimates(object)[length(object$delta), ]
    } else {
        sum$coefficients2 <-
            glarmaModelEstimates(object)[(ncol(object$X) + 1):
                                         length(object$delta), ]
    }
    class(sum) <- "summary.glarma"
    sum
}



## print method for summary method of glarma
print.summary.glarma <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 ...) {
    cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat(as.character(x$method$Resid.type), "Residuals:\n",
        sep = " ")
    if (x$df.residual >5) {
        x$deviance.resid <- setNames(quantile(x$deviance.resid, na.rm = TRUE),
                                     c("Min", "1Q", "Median", "3Q", "Max"))
    }
    xx <- zapsmall(x$deviance.resid, digits + 1L)
    print.default(xx, digits = digits, na.print = "", print.gap = 2L)
    ## print.default(format(x$five.num.summary, digits = 4), print.gap = 2,
    ##               quote = FALSE)
    if (x$methods$type == "NegBin"){
        cat("\nNegative Binomial Parameter:\n")
        printCoefmat(x$coefficients3, digits = digits, signif.legend = F, ...)
        if (x$pq > 0) {
            cat("\nGLARMA Coefficients:\n")
            printCoefmat(x$coefficients2, digits = digits,
                         signif.legend = F, ...)
        }
    } else{
        if (x$pq > 0) {
            cat("\nGLARMA Coefficients:\n")
            printCoefmat(x$coefficients2, digits = digits,
                         signif.legend = F, ...)
        }
    }
    cat("\nLinear Model Coefficients:\n")
    printCoefmat(x$coefficients1, digits = digits,
                 signif.legend = !x$tests, ...)


    cat("\n",
        apply(cbind(paste(format(c("Null", "Residual"), justify = "right"),
                          "deviance:"),
                    format(unlist(x[c("null.deviance", "deviance")]),
                           digits = max(5L, digits + 1L)), " on",
                    format(unlist(x[c("df.null", "df.residual")])),
                    " degrees of freedom\n"),
              1L, paste, collapse = " "), sep = "")
    cat("AIC:", x$aic, "\n\n", sep = " ")
    cat("Number of", switch(as.character(x$method$iter.method),
                            FS = "Fisher Scoring", NR = "Newton Raphson"),
        "iterations:", x$iter, sep = " ")
    cat("\n")
    if (x$tests == TRUE){
        cat("\nLRT and Wald Test:\n")
        cat("Alternative hypothesis: model is a GLARMA process\n")
        cat("Null hypothesis: model is a GLM with the same regression structure\n")
        printCoefmat(x$likTests, P.values = TRUE, has.Pvalue = TRUE,
                     digits = digits)
        cat("\n")
    }
}
