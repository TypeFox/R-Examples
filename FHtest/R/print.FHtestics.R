print.FHtestics <-
function (x, digits = max(options()$digits - 4, 3), ...) 
{
    saveopt <- options(digits = digits)
    on.exit(options(saveopt))
    if (!inherits(x, "FHtestics")) 
        stop("Object is not of class FHtestics")
    cat("\n")
    writeLines(x$information)
    cat("\n")
    cat(x$data.name, sep = "\n")
    cat("\n")
    otmp <- x$obs
    etmp <- x$exp
    temp <- cbind(x$n, x$diff)
    dimnames(temp) <- list(names(x$n), c("N", "O-E"))
    print(temp)
    if (substr(x$information, 2, 2) == "K") 
        cat("\nChisq= ", format(round(x$statistic, 1)), " on ", 
            length(x$n) - 1, " degrees of freedom, p-value= ", 
            format(signif(x$pvalue, digits)), "\n", sep = "")
    else cat("\nStatistic Z= ", format(round(x$statistic, 1)), 
        ", p-value= ", format(signif(x$pvalue, digits)), "\n", 
        sep = "")
    cat(x$alt.phrase, sep = "\n")
    cat("\n")
    invisible(x)
}
