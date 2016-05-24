print.noia.common <-
function (x, ...) 
{
    cat("\nPhenotype:\n")
    cat(paste("\tn=", length(x$phen), " min: ", format(min(x$phen), 
        digits = 5), " max: ", format(max(x$phen), digits = 5), 
        " mean: ", format(mean(x$phen), digits = 5), "\n"))
    cat("Genotype:\n")
    cat(paste("\tn=", nrow(x$genZ), ",", ncol(x$genZ)/3, "loci\n"))
    for (l in 1:(ncol(x$genZ)/3)) {
        cat("\t\tLocus ", l, "\t\t1: ", format(sum(x$genZ[, 3 * 
            l - 2])/(nrow(x$genZ)), digits = 3, nsmall = 3), 
            "\t2: ", format(sum(x$genZ[, 3 * l - 1])/(nrow(x$genZ)), 
                digits = 3, nsmall = 3), "\t3: ", format(sum(x$genZ[, 
                3 * l])/(nrow(x$genZ)), digits = 3, nsmall = 3), 
            "\n", sep = "")
    }
    cat("\n")
    coef <- cbind(x$E, x$variances, x$std.err, x$pvalues)
    colnames(coef) <- c("Effects", "Variances", "Std.err", "Pr(>|t|)")
    printCoefmat(coef, P.values = TRUE, signif.stars = TRUE, 
        has.Pvalue = TRUE)
    variance <- var(x$phen, na.rm = TRUE)
    cat("\nVariances\n\tTotal (phen)\t", format(variance, digits = 5), 
        "\n\tResidual\t", format(x$resvar, digits = 5), "\n\tExplained\t", 
        format(variance - x$resvar, digits = 5), "\t(", format(100 * 
            (variance - x$resvar)/variance, digits = 3), "%)", 
        "\n\tGenetic \t", format(sum(x$variances, na.rm = TRUE), 
            digits = 5), "\n", sep = "")
}
