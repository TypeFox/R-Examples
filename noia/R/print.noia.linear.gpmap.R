print.noia.linear.gpmap <-
function (x, ...) 
{
    cat("\nGenotypic values:\n")
    cat(paste("\tn=", length(x$gmap), " min: ", format(min(x$gmap), 
        digits = 5), " max: ", format(max(x$gmap), digits = 5), 
        " mean: ", format(mean(x$gmap), digits = 5), "\n"))
    cat("Genotype:\n")
    cat(paste("\t", x$nloc, " loci\n"))
    for (l in 1:x$nloc) {
        cat("\t\tLocus ", l, "\t\t1: ", format(x$genofreqloc[l, 
            1], digits = 3, nsmall = 3), "\t2: ", format(x$genofreqloc[l, 
            2], digits = 3, nsmall = 3), "\t3: ", format(x$genofreqloc[l, 
            3], digits = 3, nsmall = 3), "\n", sep = "")
    }
    cat("\n")
    coef <- cbind(x$E, x$variances)
    rownames(coef) <- names(x$E)
    colnames(coef) <- c("Effects", "Variances")
    printCoefmat(coef, P.values = FALSE, signif.stars = FALSE, 
        has.Pvalue = FALSE)
    variance <- x$V_G
    cat("\nVariances\n\tTotal (genetic)\t", format(variance, 
        digits = 5), "\n", sep = "")
}
