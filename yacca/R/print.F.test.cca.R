`print.F.test.cca` <- function(x, ...)
{
    cat("\n\tF Test for Canonical Correlations (Rao's F Approximation)\n\n")
    ctab <- cbind(x$corr, x$statistic, x$parameter, x$p.value)
    colnames(ctab) <- c("Corr","F","Num df","Den df","Pr(>F)")
    rownames(ctab) <- names(x$corr)
    printCoefmat(ctab)
    cat("\n")
}
