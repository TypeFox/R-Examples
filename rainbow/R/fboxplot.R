fboxplot = function (data, plot.type = c("functional", "bivariate"), type = c("bag", 
    "hdr"), alpha = c(0.05, 0.5), projmethod = c("PCAproj", "rapca"), factor = 1.96, na.rm = TRUE, 
    xlab = data$xname, ylab = data$yname, shadecols = gray((9:1)/10), 
    pointcol = 1, plotlegend = TRUE, legendpos = "topright", ncol=2, ...) 
{
	projmethod = match.arg(projmethod)
    op <- par(no.readonly = TRUE)
    type <- match.arg(type)
    plot.type <- match.arg(plot.type)
    y = t(data$y)
    if (na.rm) 
        y <- na.omit(y)
    if (plot.type == "functional") {
        if (type == "bag") 
            fbag(data, factor, xlab = xlab, ylab = ylab, plotlegend = plotlegend, legendpos = legendpos, ncol=ncol, projmethod = projmethod, ...)
        else fhdr(data, alpha, xlab = xlab, ylab = ylab, plotlegend = plotlegend, legendpos = legendpos, ncol=ncol, projmethod = projmethod, ...)
    }
    if (plot.type == "bivariate") {
        par(pty = "s")
        sco = PCAproj(t(data$y), center = median)$scores
        if (type == "bag") 
            bbag(data, factor, projmethod = projmethod, ...)
        else bhdr(data, alpha, shadecols = shadecols, pointcol = pointcol, projmethod = projmethod, ...)
        exit.restore <- function() {
            par(op)
        }
        on.exit(exit.restore())
    }
}

