plot.fscores <-
function (x, bw = "nrd0", adjust = 2, kernel = "gaussian", 
            include.items = FALSE, tol = 0.2, xlab = "Ability", ylab = "Density", 
            main = "Kernel Density Estimation for Ability Estimates", pch = 16, cex = 1.5, ...) {
    if (!inherits(x, "fscores"))
        stop("Use only with 'fscores' objects.\n")
    if (x$resp.pats)
        stop("the plot() method requires 'resp.patterns = NULL'.\n")
    ablts <- rep(x$score.dat$z1, x$score.dat$Obs)
    dens <- density(ablts, bw = bw, adjust = adjust, kernel = kernel)
    if (include.items) {
        if (is.null(x$coef)) {
            warning("'include.items' works only for dichotomous data.\n")
            plot(dens, xlab = xlab, ylab = ylab, main = main, ...)
        } else {
            if (!any((cnams <- colnames(x$coef)) == "Dffclt"))
                warning("when 'include.items = TRUE' you should fit the model under the IRT parameterization.\n")
            betas <- if (cnams[1] == "Gussng" || cnams[1] == "c.i") x$coef[, 2] else x$coef[, 1]
            plot(dens, xlim = range(dens$x, betas), xlab = xlab, ylab = ylab, main = main, ...)
            Beta <- round(betas / tol) * tol
            stripchart(Beta, at = 0, method = "stack", cex = cex, pch = pch, add = TRUE)
        }
    } else {
        plot(dens, xlab = xlab, ylab = ylab, main = main, ...)
    }
    invisible()
}
