xyplot.z <-
function (xx, yy, zz, p = 0.5, sfact = 2.5, zmin = NA, zmax = NA, 
    log = NULL, xlim = NULL, ylim = NULL, xlab = deparse(substitute(xx)), 
    ylab = deparse(substitute(yy)), zlab = deparse(substitute(zz)), 
    main = "", col = 1, iflgnd = FALSE, title = deparse(substitute(zz)), 
    cex.lgnd = 0.8, ifparams = FALSE, cex.params = 0.8, ...) 
{
    frame()
    temp.z <- remove.na(cbind(xx, yy, zz))
    x <- temp.z$x[1:temp.z$n, 1]
    y <- temp.z$x[1:temp.z$n, 2]
    z <- temp.z$x[1:temp.z$n, 3]
    nz <- temp.z$n
    if (main == "") 
        if (zlab == "") 
            banner <- ""
        else banner <- paste("Proportional Symbol Plot for", 
            zlab)
    else banner <- main
    z.min <- min(z)
    z.max <- max(z)
    zrange <- c(zmin, zmax)
    rgz <- syms(z, zrange, p = p)
    if (is.null(log)) 
        log <- ""
    plot(x, y, type = "n", xlab = xlab, ylab = ylab, xlim = xlim, 
        ylim = ylim, log = log, main = banner, ...)
    symbols(x, y, circles = rgz, inches = sfact * 0.05, fg = col, 
        add = TRUE, ...)
    if (iflgnd) {
        if (!is.na(zmax)) 
            z[z > zmax] <- zmax
        if (!is.na(zmin)) 
            z[z < zmin] <- zmin
        zval <- quantile(z, prob = c(1, 0.75, 0.5, 0.25, 0))
        rgz <- syms(zval, zrange, p = p)
        zval <- signif(zval, 3)
        legend(locator(1), pch = rep.int(1, 5), pt.cex = rgz[1:5] * 
            sfact*1.25, col = rep.int(col, 5), paste(" ", 
            zval[1:5]), cex = cex.lgnd, title = title)
    }
    if (ifparams) 
        text(locator(1), paste("p =", signif(p, 3), "& sfact =", 
            sfact, "\nz.max =", signif(z.max, 3), "; zmax =", 
            zmax, "\nz.min =", signif(z.min, 3), "; zmin =", 
            zmin), adj = c(0, 1), cex = cex.params, ...)
    invisible()
}
