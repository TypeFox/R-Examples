`plot.ftsf` <- function (x, plot.type = c("function", "components", "variance"), 
    components, xlab1 = fit$y$xname, ylab1 = "Basis function", 
    xlab2 = "Time", ylab2 = "Coefficient", mean.lab = "Mean", 
    level.lab = "Level", main.title = "Main effects", interaction.title = "Interaction", 
    vcol = 1:3, shadecols = 7, fcol = 4, basiscol = 1, coeffcol = 1, 
    outlier.col = 2, outlier.pch = 19, outlier.cex = 0.5, ...) 
{
    plot.type <- match.arg(plot.type)
    if (plot.type == "function") 
        plot(x$mean, ...)
    else if (plot.type == "variance") {
        xx <- x[[1]]$x
        xname <- x[[1]]$xname
        plot(x[[1]]$x, rowMeans(x$var$model), type = "l", xlab = xname, 
            ylab = "Variance", col = vcol[1])
        abline(0, 0, lty = 2)
        lines(xx, x$var$mean, col = vcol[2])
        lines(xx, x$var$error, col = vcol[3])
    }
    else {
        fit <- x$model
        oldpar <- par(no.readonly = TRUE)
        mean <- is.element("mean", colnames(fit$basis))
        level <- is.element("level", colnames(fit$basis))
        m <- mean + level
        order <- ncol(fit$basis) - m
        if (missing(components)) 
            components <- order
        n <- components + (mean | level)
        par(mfcol = c(2, n))
        if (mean) 
            plot(fit$y$x, fit$basis[, "mean"], type = "l", lty = 1, 
                xlab = xlab1, ylab = mean.lab, main = main.title, 
                col = basiscol, ...)
        if (level) 
            plot.ts(fit$coeff[, "level"], xlab = xlab2, ylab = level.lab, 
                    main = ifelse(mean, "", main.title), col = coeffcol, 
                    ...)
        if (m == 1) 
            plot(0, 0, type = "n", xaxt = "n", yaxt = "n", bty = "n", 
                 xlab = "", ylab = "")
        if (components > 0) {
            for (i in 1:components) {
                 yl1 <- ifelse(components > 1, paste(ylab1, i), ylab1)
                 yl2 <- ifelse(components > 1, paste(ylab2, i), ylab2)
                 plot(fit$y$x, fit$basis[, m + i], type = "l", 
                      lty = 1, xlab = xlab1, ylab = yl1, col = basiscol, 
                      ...)
                 if (i == 1) 
                     title(interaction.title)
                 plot(x$coeff[[i + m]], type = "n", xlab = xlab2, 
                      ylab = yl2, main = "", shadecols = shadecols, 
                      fcol = fcol, col = coeffcol, ...)
                 lines(x$coeff[[i + m]]$x, col = "gray")
                 junk <- x$coeff[[i + m]]$x
                 junk[fit$weights < 1e-06] <- NA
                 lines(junk, col = coeffcol)
                 if (!is.null(fit$weights)) 
                     points(time(fit$coeff)[fit$weights < 0.1], 
                     fit$coeff[fit$weights < 0.1, i + m], pch = outlier.pch, 
                    col = outlier.col, cex = outlier.cex)
            }
        }
        par(oldpar)
    }
}
