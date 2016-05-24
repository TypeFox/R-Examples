`print.mantel.pertables` <-
function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("Mantel statistic based on", x$mantel$mantel.raw$method, 
        "\n\n")
    cat("Mantel statistic r: ")
    cat(formatC(-x$mantel$mantel.raw$statistic, digits = digits), 
        "\n")
    nperm <- x$mantel$mantel.raw$permutations
    if (nperm) {
        cat("      Significance:", format.pval(x$mantel$mantel.raw$signif, 
            eps = 1/nperm), "\n\n")
    }
    n <- x$simulation$mantel.quant
    cat("Confidence intervals of statistic r and p-values under different taxonomic scenarios", 
        "\n\n")
    print(n)
    cat("\n")
    cat("Significance of Mantel statistic r under different taxonomic scenarios", 
        "\n")
    cat("Significance:", format.pval(x$mantel$ptax), "\n\n")
}

