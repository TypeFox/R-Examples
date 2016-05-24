"print.permtstBurSt" <-
function (x, digits = 4, ...) 
{
    cat("Call:\n")
    print(x$call)
    cat("\nOriginal value:", x$original.score, "  Number of observations:", 
        x$stats["nobs"], "\n")
    cat("Number of random permutations:", x$stats["trials"], 
        "  Alternative:", x$alternative, "  p-value:", 
	   round(x$stats["p.value"], digits), "\n")
    invisible(x)
}
