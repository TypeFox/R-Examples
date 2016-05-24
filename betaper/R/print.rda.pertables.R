`print.rda.pertables` <-
function (x, ...) 
{
    n <- x$simulation$rda.quant[, -3]
    cat("\n")
    cat("Confidence intervals of R-squared and  pseudo-F values for RDA under different taxonomic scenarios", 
        "\n\n")
    print(n)
}

