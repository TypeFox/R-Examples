gx.rma <-
function (xx1, xx2, x1lab = NULL, x2lab = NULL, log = FALSE, 
    if.plot = FALSE, if.rma = FALSE, if.coeffs = FALSE, ...) 
{
    if (!if.plot) 
        if.rma <- FALSE
    if (is.matrix(xx1)) {
        xx1.name <- dimnames(xx1)[[2]][1]
        xx2.name <- dimnames(xx1)[[2]][2]
        x1 <- xx1[, 1]
        x2 <- xx1[, 2]
        xx1 <- x1
        xx2 <- x2
    }
    else {
        xx1.name <- deparse(substitute(xx1))
        xx2.name <- deparse(substitute(xx2))
    }
    if (length(xx1) != length(xx2)) 
        stop("Input vectors must be of equal length\n")
    temp.x <- remove.na(cbind(xx1, xx2))
    x1 <- temp.x$x[1:temp.x$n, 1]
    x2 <- temp.x$x[1:temp.x$n, 2]
    if (is.null(x1lab)) 
        x1lab <- xx1.name
    if (is.null(x2lab)) 
        x2lab <- xx2.name
    cat("\n Reduced Major Axis for", x1lab, "and", x2lab)
    if (if.plot) {
        x.min <- min(min(x1), min(x2))
        x.max <- max(max(x1), max(x2))
        if (log) 
            log.plot <- "xy"
        else log.plot <- ""
        if (is.null(x1lab)) 
            x1lab <- xx1.name
        if (is.null(x2lab)) 
            x2lab <- xx2.name
        plot(x1, x2, log = log.plot, xlim = c(x.min, x.max), 
            ylim = c(x.min, x.max), xlab = x1lab, ylab = x2lab, 
            ...)
        abline(0, 1, lty = 2)
    }
    if (log) {
        if (min(min(x1), min(x2)) <= 0) 
            stop("Vector(s) contain one or more <= 0 values\n")
        x1 <- log10(x1)
        x2 <- log10(x2)
        cat("\n Data have been Log10 transformed")
    }
    xlen <- temp.x$n
    x <- cbind(x1, x2)
    xbar <- cbind(mean(x1), mean(x2))
    xvar <- cbind(var(x1), var(x2))
    xsdv <- sqrt(xvar)
    slope <- xsdv[2]/xsdv[1]
    r <- cov(x1, x2)/(xsdv[1] * xsdv[2])
    if (r < 0) 
        slope <- slope * (-1)
    incpt <- xbar[2] - xbar[1] * slope
    temp <- (1 - r * r)/xlen
    seslp <- slope * sqrt(temp)
    seint <- xsdv[2] * sqrt(temp * (1 + (xbar[1] * xbar[1])/xvar[1]))
    temp <- qt(0.975, xlen - 1)
    cislp <- seslp * temp
    slpll <- slope - cislp
    slpul <- slope + cislp
    ciint <- seint * temp
    intll <- incpt - ciint
    intul <- incpt + ciint
    cat("\n\n Means =\t", signif(xbar[1], 4), "\t\t", signif(xbar[2], 4), 
        "\n SDs =\t\t", signif(xsdv[1], 4), "\t\t", signif(xsdv[2], 4), 
        "\n\n Corr =\t\t", round(r, 4), "\n N =\t\t", xlen, 
        "\n\t\t\t\t   SE\t\t\t95% CLs", "\n Slope =\t", signif(slope, 4),
        "\t", signif(seslp, 4), "\t  ", signif(slpll, 4), "<->",
        signif(slpul, 4), "\n Intercept =\t", signif(incpt, 4), "\t", 
        signif(seint, 6), "\t  ", signif(intll, 4), "<->", signif(intul, 4), 
        "\n\n")
    H0text <- "Reject"
    if ((abs(slope - 1) <= cislp) & (incpt <= ciint)) H0text <- "Accept"
    cat(paste(" ", H0text, " hypothesis that RMA is (0,1)\n\n", sep = ""))
    if (if.rma) 
        abline(incpt, slope, lty = 3, col = 2)
    if (if.coeffs) {
        text <- paste("Reduced Major Axis - Orthogonal Regression\n",
            "N = ", xlen, ", Coefficients:\na0 = ", signif(incpt, 3),
            "  95% CLs: ", signif(intll, 4), "<->", signif(intul, 4),
            "\na1 = ", signif(slope, 4), "  95% CLs: ",
            signif(slpll, 4), "<->", signif(slpul, 4), "\n", H0text,
            " hypothesis that RMA is (0,1)", sep = "")
        text(locator(1), text, adj = c(0, 1), cex = 0.8)
    }
    invisible(list(n = xlen, mean = xbar, sd = xsdv, corr = r, 
        a0 = incpt, a1 = slope, sea0 = seint, sea1 = seslp))
}
