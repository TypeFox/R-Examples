gx.triples.aov <-
function (x, xname = deparse(substitute(x)), log = FALSE, table = FALSE) 
{
    n <- length(x)
    if (n != length(na.omit(x))) 
        stop("Missing data, NAs, not permitted")
    ntrip <- round(n/3)
    if (n != ntrip * 3) 
        stop("Incomplete set of triplicates")
    if (log) {
        data.name <- xname
        x <- log10(x)
        cat("\n Data for", data.name, "have been log10 tansformed\n")
    }
    xx <- matrix(x, nrow = ntrip, ncol = 3, byrow = TRUE)
    if (table) {
        cat("\n Input triplicate data for", xname, "\n")
        for (i in 1:ntrip) cat(" ", xx[i, 1:3], "\n")
    }
    ba.ss <- sum((xx[, 1] - xx[, 2])^2)/2
    ba.df <- ws.df <- ntrip
    ba.ms <- ba.ss/ba.df
    bs.df <- ntrip - 1
    bs.ss <- var(rowMeans(xx)) * 3 * bs.df
    bs.ms <- bs.ss/bs.df
    t.df <- n - 1
    t.ss <- var(x) * t.df
    ws.ss <- t.ss - (bs.ss + ba.ss)
    ws.ms <- ws.ss/ws.df
    vcomp.pct <- vcomp <- numeric(3)
    vcomp[1] <- ba.ms
    vcomp[2] <- (ws.ms - ba.ms)/1.3333
    vcomp[3] <- (bs.ms - (ba.ms + 1.6667 * vcomp[2]))/3
    F.prob <- numeric(2)
    F1 <- ws.ms/ba.ms
    F.prob[1] <- pf(F1, ws.df, ba.df)
    swsms <- vcomp[1] + vcomp[2] * 1.6667
    swsdf <- swsms^2/((vcomp[1]^2)/ba.df + ((vcomp[2] * 1.6667)^2)/ws.df)
    F2 <- bs.ms/swsms
    F.prob[2] <- pf(F2, bs.df, swsdf)
    p.char <- character(2)
    for (k in 1:2) {
        if (F.prob[k] > 0.9999) 
            F.prob[k] <- 0.9999
        p.char[k] <- " NS"
        if (F.prob[k] > 0.95) 
            p.char[k] <- "  *"
        if (F.prob[k] > 0.99) 
            p.char[k] <- " **"
        if (F.prob[k] > 0.999) 
            p.char[k] <- "***"
    }
    for (k in 1:3) if (vcomp[k] < 0) 
        vcomp[k] = 0
    totvar <- sum(vcomp)
    vcomp.pct <- 100 * vcomp/totvar
    cat("\n ANOVA ('Bainbridge' staggered design) for:", xname, "\n")
    cat("\n Level   Variation\tSums of     DF    Mean\t\t Synthesized (1)     F       Prob (2)", 
        "\n\t  between\tSquares          Square\t\t   MS       DF     Ratio     Value\n")
    cat("\n   3     Sites\t\t", signif(bs.ss, 4), "\t    ", bs.df, 
        "\t", signif(bs.ms, 5), "\t", signif(swsms, 5), "  ", 
        signif(swsdf, 3), "  ", signif(F2, 4), "   ", signif(F.prob[2], 
            4), p.char[2], "\n   2     Samples @ 3\t", signif(ws.ss, 
            4), "\t    ", ws.df, "\t", signif(ws.ms, 5), "\t\t\t  ", 
        signif(F1, 4), "   ", signif(F.prob[1], 4), p.char[1], 
        "\n   1     Analyses\t", signif(ba.ss, 4), "\t    ", 
        ba.df, "\t", signif(ba.ms, 5), "\n\n Totals\t\t\t", signif(t.ss, 
            4), "\t   ", t.df, "\n", "\n Notes: (1) Synthesized MS and DF calculated after Satterthwaite (1946)", 
        "\n        (2) NS = Not Significant, F <= 0.95\n\t     * = F > 0.95 & <= 0.99", 
        "\n\t    ** = F > 0.99 & <= 0.999\n\t   *** = F > 0.999\n")
    cat("\n Level   Variation\tSums of     DF    Mean    Unit    Variance       %age", 
        "\n\t  between\tSquares          Square   Size    Component\n")
    cat("\n   3     Sites\t\t", signif(bs.ss, 4), "\t    ", bs.df, 
        "\t", signif(bs.ms, 5), "  ", ntrip, "\t", signif(vcomp[3], 
            5), "\t", signif(vcomp.pct[3], 3), "\n   2     Samples @ 3\t", 
        signif(ws.ss, 4), "\t    ", ws.df, "\t", signif(ws.ms, 
            5), "  ", ntrip * 2, "\t  ", signif(vcomp[2], 5), 
        "\t", signif(vcomp.pct[2], 3), "\n   1     Analyses\t", 
        signif(ba.ss, 4), "\t    ", ba.df, "\t", signif(ba.ms, 
            5), "  ", ntrip * 3, "\t  ", signif(vcomp[1], 5), 
        "\t", signif(vcomp.pct[1], 3), "\n\n Totals\t\t\t", signif(t.ss, 
            4), "\t   ", t.df, "\t\t\t  ", signif(totvar, 5), 
        "\n")
    v <- vcomp[3]/(vcomp[2] + vcomp[1])
    cat("\n Empirical variance ratio, V, a measure of sampling efficiency =", 
        round(v, 2), "\n (Note: V>>1 is desirable)\n\n")
    invisible()
}
