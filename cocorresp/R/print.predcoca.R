`print.predcoca` <- function(x, digits = max(3, getOption("digits") - 3),
                             ...) {
    cat("\nPredictive Co-Correspondence Analysis\n\n")
    writeLines(strwrap(pasteCall(x$call)))
    if (!is.null(x$lambda)) {
        cat("\nEigenvalues:\n")
        print(round(x$lambda, digits), ..., print.gap = 2)
    } else {
        pcentX <- (x$varianceExp$Xblock / x$totalVar$Xblock) * 100
        Xvar.mat <- rbind(pcentX, cumsum(pcentX))
        pcentY <- (x$varianceExp$Yblock / x$totalVar$Yblock) * 100
        Yvar.mat <- rbind(pcentY, cumsum(pcentY))
        rownames(Xvar.mat) <- rownames(Yvar.mat) <- c("Individual:", "Cumulative:")
        cat("\nPercentage Variance Explained:\n")
        cat("\nY-block: variance explained in", x$nam.dat$namY,
            "(response) \n", sep = " ")
        print(round(Yvar.mat, digits), ..., print.gap = 2)
        cat("\nX-block: variance explained in", x$nam.dat$namX,
            "(predictor) \n", sep = " ")
        print(round(Xvar.mat, digits), ..., print.gap = 2)
    }
    cat("\n")
    invisible(x)
}
