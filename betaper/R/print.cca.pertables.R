`print.cca.pertables` <-
function (x, ...) 
{
    n <- x$simulation$cca.quant[, -3]
    cat("\n")
    cat("Confidence intervals of R-squared and  pseudo-F values for CCA under different taxonomic scenarios", 
        "\n\n")
    print(n)
}

