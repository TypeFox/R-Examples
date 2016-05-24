# data ellipses for xyplot
# M. Friendly, 7/7/2010 9:24AM
# - replaced panel.polygon with panel.lines

# calculate one data ellipse; modified from car to do just the calculation

conf.ellipse <- function(x, y, level = 0.68, segments = 50L, robust = FALSE)
{
    if (length(x) != length(y))
        stop("x and y must be vectors of the same length")
    xy <- na.omit(cbind(x, y))
    if (robust)
    {
        v <- MASS::cov.trob(xy)
        shape <- v$cov
        center <- v$center
    }
    else
    {
        shape <- var(xy)
        center <- colMeans(xy)
    }
    radius <- sqrt(2 * qf(level, df1 = 2, df2 = length(x) - 1))
    angles <- seq(0, 2 * pi, length.out = segments + 1)
    unit.circle <- cbind(cos(angles), sin(angles))
    ell <- t(center + radius * t(unit.circle %*% chol(shape)))
    list(ellipse = ell, center = center)
}

panel.ellipse <-
    function(x, y, groups = NULL,
             level = 0.68, segments = 50, robust = FALSE,
             center.pch = 3, center.cex = 2, 
             ..., type, pch, cex)
{
    x <- as.numeric(x)
    y <- as.numeric(y)
    if (!is.null(groups))
        panel.superpose(x = x, y = y, groups = groups,
                        level = level, segments = segments, robust = robust,
                        center.pch = center.pch, center.cex = center.cex,
                        panel.groups = panel.ellipse, ...)
    else
    {
        ell <- conf.ellipse(x = x, y = y, level = level, segments = segments, robust = robust)
        panel.xyplot(ell$ellipse[,1], ell$ellipse[,2], ..., type = "l")
        if (!is.null(center.pch))
            panel.xyplot(ell$center[1], ell$center[2], ...,
                         pch = center.pch, cex = center.cex)
    }
}


if (FALSE)
{

    old.panel.ellipse <-
        function(x, y, groups=NULL, level=0.68, segments=51, robust=FALSE,
                 col=NA, col.line = if (is.null(groups)) plot.line$col else superpose.line$col,
                 lwd = if (is.null(groups)) plot.line$lwd else superpose.line$lwd, 
                 lty = if (is.null(groups)) plot.line$lty else superpose.line$lty,
                 center.pch=3, center.cex=2, 
                 ...)
        {
            x <- as.numeric(x)
            y <- as.numeric(y)
            plot.line <- trellis.par.get("plot.line")
            superpose.line <- trellis.par.get("superpose.line")

            if (!is.na(col))
            {
                if (missing(col.line)) col.line <- col
            }
            
            groups <- as.factor(if(is.null(groups)) rep(1, length(x)) else as.character(groups))
            n.groups <- length(levels(groups))
            col <- rep(col.line, length.out=n.groups)
            lty <- rep(lty, length.out=n.groups)
            lwd <- rep(lwd, length.out=n.groups)
            for (i in 1:n.groups)
            {
                subs <- groups == levels(groups)[i]
                XY <- na.omit(data.frame(x=x[subs], y=y[subs]))
                ell <- .ellipse(XY, level = level,
                                segments = segments,
                                robust = robust)
                ## panel.polygon(ell[,1], ell[,2], border=col[i], col="transparent", lty=lty[i], lwd=lwd[i], ...)
                panel.lines(ell[,1], ell[,2], col = col[i], lty = lty[i], lwd = lwd[i], ...)
                if (!is.null(center.pch)) {
                    center = colMeans(XY)
                    panel.points(center[1], center[2], col=col[i],
                                 pch = center.pch, cex = center.cex)
                }
            }
        }

}




