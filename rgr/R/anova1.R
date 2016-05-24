anova1 <-
function (x1, x2, xname = deparse(substitute(x1)), log = FALSE) 
{
    if (length(x1) != length(x2)) {
        cat("\nLengths of vectors not the same\n")
        return()
    }
    temp.x <- remove.na(cbind(x1, x2))
    a1 <- temp.x$x[1:temp.x$n, 1]
    a2 <- temp.x$x[1:temp.x$n, 2]
    cat("\n Combined Sampling and Analytical, or Analytical Variability, Study,", 
        "\n Utilizes Field Sampling or Laboratory Duplicates.  In ANOVA Tables,", 
        "\n the variability:", "\n   Between would be between sampling sites or analysed samples, and", 
        "\n   Within would be at sampling sites or due to duplicate analyses")
    if (log) {
        if (min(min(a1), min(a2)) <= 0) {
            cat("\nVector(s) contain one or more <= 0 values\n")
            return()
        }
        a1 <- log10(a1)
        a2 <- log10(a2)
        cat("\n\n Data have been Log10 transformed for the ANOVA")
    }
    alen <- length(a1)
    a <- c(a1, a2)
    adiff <- a1 - a2
    sams <- sum(adiff * adiff)/(2 * alen)
    tms <- var(a)
    tdf <- alen * 2 - 1
    tss <- tms * tdf
    rss <- 2 * var(adiff)
    rdf <- alen - 1
    rms <- rss/rdf
    wss <- (sams * alen) - rss
    wdf <- 1
    bss <- tss - (wss + rss)
    bdf <- alen - 1
    bms <- bss/bdf
    fval1 <- bms/rms
    prob1 <- 1 - pf(fval1, bdf, 1)
    fval2 <- wss/rms
    prob2 <- 1 - pf(fval2, 1, rdf)
    cat("\n\n Two-Way Random Effects Model for", xname, "\n Source\t\t  SS\t\tdf\t  MS\t\t  F\t Prob", 
        "\n Between\t", format(signif(bss, 5)), "\t", bdf, "\t", 
        format(signif(bms, 5)), "\t\t", format(round(fval1, 2)), 
        "\t", format(round(prob1, 4)), "\n Within\t\t", format(signif(wss, 
            5)), "\t\t", wdf, "\t", format(signif(wss, 5)), "\t\t", 
        format(round(fval2, 2)), "\t", format(round(prob2, 4)), 
        "\n Residual\t", format(signif(rss, 5)), "\t\t", rdf, 
        "\t", format(signif(rms, 5)), "\n Total\t\t", format(signif(tss, 
            5)), "\t", tdf, "\t", format(signif(tms, 5)))
    wss <- sams * alen
    bss <- tss - wss
    bdf <- alen - 1
    bms <- bss/bdf
    fval1 <- round(bms/sams, 2)
    prob1 <- 1 - pf(fval1, bdf, alen)
    cat("\n\n One-Way Random Effects Model for", xname, "\n Source\t\t  SS\t\tdf\t  MS\t\t  F\t Prob", 
        "\n Between\t", format(signif(bss, 5)), "\t", bdf, "\t", 
        format(signif(bms, 5)), "\t\t", format(round(fval1, 2)), 
        "\t", format(round(prob1, 4)), "\n Within\t\t", format(signif(wss, 
            5)), "\t\t", alen, "\t", format(signif(sams, 5)), 
        "\n Total\t\t", format(signif(tss, 5)), "\t", tdf, "\t", 
        format(signif(tms, 5)))
    bvar <- (bms - sams)/2
    tvar <- bvar + sams
    bpct <- (100 * bvar)/tvar
    wpct <- (100 * sams)/tvar
    cat("\n\n Source\t\t  MS\t\tVar Comp\t %age", "\n Between\t", 
        format(signif(bms, 5)), "\t\t", format(signif(bvar, 5)), 
        "\t\t", format(round(bpct, 1)), "\n Within\t\t", format(signif(sams, 
            5)), "\t\t", format(signif(sams, 5)), "\t\t", format(round(wpct, 
            1)), "\n\t\t\t\t", format(signif(tvar, 5)))
    v <- round(bvar/sams, 2)
    vm <- round((2 * bvar)/sams, 2)
    xmean <- mean(a)
    rsd <- (100 * sqrt(sams))/xmean
    cat("\n\n Summary Statistics for", xname, "\n Grand Mean =\t", 
        format(signif(xmean, 5)), "\t\t Variance =\t", format(signif(tms, 
            5)), "\n 'Error' S^2 =\t", format(signif(sams, 5)), 
        "\t\t Std. Dev. =\t", format(signif(sqrt(sams), 5)), 
        "\n 'Error' RSD% =\t", format(round(rsd, 1)), "\n Miesch's V = \t", 
        v, "\t\t Vm = \t\t", vm, "\n\n")
    invisible()
}
