"print.summary.crossval" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    ##calculations/manipulations
    names(x$CVfit) <- paste("COCA", 1:x$n.axes, sep = " ")
    pcentX <- (x$varianceExp$Xblock / x$totalVar$Xblock) * 100
    Xvar.mat <- rbind(pcentX, cumsum(pcentX))
    pcentY <- (x$varianceExp$Yblock / x$totalVar$Yblock) * 100
    Yvar.mat <- rbind(pcentY, cumsum(pcentY))
    rownames(Xvar.mat) <- rownames(Yvar.mat) <- c("Individual:", "Cumulative:")
    colnames(Xvar.mat) <- colnames(Yvar.mat) <- names(x$CVfit)
    ##printing
    cat("\nCross-validation for Predictive Co-Correspondence Analysis\n\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat(sprintf("\nCross-validatory %%fit of %s to %s:\n\n", x$nam.dat$namX,
                x$nam.dat$namY))
    print(round(x$CVfit, digits), print.gap = 2)
    cat(sprintf("\nTotal Variance in %s: %f\n", x$nam.dat$namY,
                format(x$totalVar$Yblock, digits)))
    cat(sprintf("\nTotal Variance in %s: %f\n", x$nam.dat$namX,
                format(x$totalVar$Xblock, digits)))
    cat("\nPercentage Variance Explained:\n")
    cat("\nY-block: variance explained in", x$nam.dat$namY,
        "(response) \n", sep = " ")
    print(round(Yvar.mat, digits), ..., print.gap = 2)
    cat("\nX-block: variance explained in", x$nam.dat$namX,
        "(predictor) \n", sep = " ")
    print(round(Xvar.mat, digits), ..., print.gap = 2)
    invisible(x)
  }

