"print.permutest.coca" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    #Alter this to only display the axes and thier p-values
    # The summary method will display the full matrix
    ptest.stats <- rbind(x$permstat, x$inertia,
                         x$fitax, x$pcent.fit, x$pval)
    rownames(ptest.stats) <- c("Stat.", "Inertia",
                               "Fit", "% fit", "P-value")
    colnames(ptest.stats) <- paste("COCA", 1:x$n.axes, sep = " ")
    cat("\nPermutation test for predictive co-correspondence analysis:\n\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat("\nPermutation test results:\n\n")
    printCoefmat(t(ptest.stats), digits = digits, na.print = "")
    cat("\n")
    invisible(x)
  }

