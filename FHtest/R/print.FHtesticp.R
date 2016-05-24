print.FHtesticp <-
function (x, digits = max(options()$digits - 4, 3), ...) 
{
    saveopt <- options(digits = digits)
    on.exit(options(saveopt))
    if (!inherits(x, "FHtesticp")) 
        stop("Object is not of class FHtesticp")
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
    if (length(grep("exact", x$information)) == 1) 
        cat("\nStatistic= ", format(round(x$statistic, 1)), ", p-value= ", 
            format(signif(x$pvalue, digits)), "\n", sep = "")
    else if (substr(x$information, 2, 2) == "K") 
        cat("\nChisq= ", format(round(x$statistic, 1)), " on ", 
            length(x$n) - 1, " degrees of freedom, p-value= ", 
            format(signif(x$pvalue, digits)), "\n", sep = "")
    else cat("\nStatistic Z= ", format(round(x$statistic, 1)), 
        ", p-value= ", format(signif(x$pvalue, digits)), "\n", 
        sep = "")
    if (length(grep("Monte", x$information)) == 1) 
        cat(format(100 * attr(x$p.conf.int, "conf.level")), "percent confidence interval on p-value: ", 
            format(c(x$p.conf.int[1], x$p.conf.int[2])), "\n")
    cat(x$alt.phrase, sep = "\n")
    cat("\n")
    invisible(x)
}
