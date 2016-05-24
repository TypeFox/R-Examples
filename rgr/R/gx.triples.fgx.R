gx.triples.fgx <-
function (x, RepStat, xname = deparse(substitute(x)), log = FALSE) 
{
    temp.x <- remove.na(cbind(x, RepStat))
    xx <- temp.x$x[1:temp.x$n, 1]
    RS <- temp.x$x[1:temp.x$n, 2]
    if (log) {
        data.name <- deparse(substitute(x))
        xx <- log10(xx)
        cat(paste("\n Data for", data.name, "have been log10 tansformed"))
    }
    cat("\n F-tests to check on triples representivity for:\n", 
        xname)
    var0 <- var(xx[RS == 0])
    df0 <- length(xx[RS == 0]) - 1
    var1 <- var(xx[RS == 1])
    df1 <- length(xx[RS == 1]) - 1
    var2 <- var(xx[RS == 2])
    df2 <- length(xx[RS == 2]) - 1
    F.prob <- numeric(2)
    F0 <- var0/var1
    if (F0 > 1) {
        F.prob[1] <- pf(F0, df0, df1)
        df11 <- df0
        df12 <- df1
        note <- " (Regional > triples)"
    }
    else {
        F0 <- 1/F0
        F.prob[1] <- pf(F0, df1, df0)
        df11 <- df1
        df12 <- df0
        note <- " (Regional < triples)"
    }
    F12 <- var1/var2
    if (F12 > 1) {
        F.prob[2] <- pf(F12, df1, df2)
        df21 <- df1
        df22 <- df2
    }
    else {
        F12 <- 1/F12
        F.prob[2] <- pf(F12, df2, df1)
        df21 <- df2
        df22 <- df1
    }
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
    cat("\n\n F(regional) =", signif(F0, 4), "\tDoFs =", df11, 
        "and", df12, "\tProb =", signif(F.prob[1], 3), p.char[1], 
        note)
    cat("\n F(1 vs. 2)  =", signif(F12, 4), "\tDoFs =", df21, 
        "and", df22, "\t\tProb =", signif(F.prob[2], 3), p.char[2], 
        "\n\n")
    invisible()
}
