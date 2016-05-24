"print.summary.permutest.coca" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    #needs an axes argument for displaying the chi-square residuals
    ptest.stats <- rbind(x$permstat, x$inertia,
                         x$fitax, x$pcent.fit, x$pval)
    rownames(ptest.stats) <- c("Stat.", "Inertia",
                               "Fit", "% fit", "P-value")
    colnames(ptest.stats) <- paste("COCA", 1:x$n.axes, sep = " ")
    cat("\nPermutation test for predictive co-correspondence analysis:\n\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat("\nTotal Inertia: ", x$total.inertia, "\n")
    cat("\nPermutation test results:\n")
    printCoefmat(t(ptest.stats), digits = digits, na.print = "")
    #cat("\nChi-square residual matrix for Y1:\n")
    #print(x$Ychi1)
    #cat("\n\nChi-square residual matrix for Y2:\n")
    #print(x$Ychi2)
    invisible(x)
  }

