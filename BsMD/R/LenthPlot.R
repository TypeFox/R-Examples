LenthPlot <-
function (obj, alpha = 0.05, plt = TRUE, limits = TRUE, xlab = "factors", 
    ylab = "effects", faclab = NULL, cex.fac = par("cex.lab"), 
    cex.axis = par("cex.axis"), adj = 1, ...) 
{
    if (inherits(obj, "lm")) {
        i <- pmatch("(Intercept)", names(coef(obj)))
        if (!is.na(i)) 
            obj <- 2 * coef(obj)[-pmatch("(Intercept)", names(coef(obj)))]
    }
    b <- obj
    if (!is.null(faclab)) {
        if (!is.list(faclab)) 
            stop("* Argument 'faclab' has to be NULL or a list with 'idx' and 'lab' elements")
        names(b) <- rep("", length(b))
        names(b)[faclab$idx] <- faclab$lab
    }
    m <- length(b)
    d <- m/3
    s0 <- 1.5 * median(abs(b))
    cj <- as.numeric(b[abs(b) < 2.5 * s0])
    PSE <- 1.5 * median(abs(cj))
    ME <- qt(1 - alpha/2, d) * PSE
    gamma <- (1 + (1 - alpha)^(1/m))/2
    SME <- qt(gamma, d) * PSE
    if (plt) {
        n <- length(b)
        x <- seq(n)
        ylim <- range(c(b, 1.2 * c(ME, -ME)))
        plot(x, b, xlim = c(1, n + 1), ylim = ylim, type = "n", 
            xlab = xlab, ylab = ylab, frame = FALSE, axes = FALSE, 
            ...)
        idx <- x[names(b) != ""]
        text(x[idx], rep(par("usr")[3], length(idx)), labels = names(b)[idx], 
            cex = cex.fac, xpd = NA)
        axis(2, cex.axis = cex.axis)
        for (i in seq(along = x)) segments(x[i], 0, x[i], b[i], 
            lwd = 3, col = 1, lty = 1)
        abline(h = 0, lty = 4, xpd = FALSE)
        if (limits) {
            abline(h = ME * c(1, -1), xpd = FALSE, lty = 2, col = grey(0.2))
            text(adj * (n + 1) * c(1, 1), (ME + strheight("M", 
                cex = cex.axis)) * c(1, -1), labels = "ME", cex = 0.9 * 
                cex.axis, xpd = FALSE)
            abline(h = SME * c(1, -1), xpd = FALSE, lty = 3, 
                col = grey(0.2))
            text(adj * (n + 1) * c(1, 1), (SME + strheight("M", 
                cex = cex.axis)) * c(1, -1), labels = "SME", 
                cex = 0.9 * cex.axis, xpd = FALSE)
        }
    }
    return(c(alpha = alpha, PSE = PSE, ME = ME, SME = SME))
}
