print.summary.StratSel <-
function(x, ...){
    print(x$call)
    cat("\n")
    printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE,  ...)
    writeLines("---")
	writeLines(paste("Number of Observations:", x$df, "-- Number of Iterations:", x$nits, sep=" "))
	writeLines(paste("AIC:", x$AIC, "-- AIC (small sample correction):", x$AIC.c, sep=" "))
	writeLines(paste("Log-Likelihood:", round(x$ll,3), "-- Convergence code:", x$conv, sep=" "))
    cat("\n")
 }
