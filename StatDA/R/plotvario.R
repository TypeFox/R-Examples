"plotvario" <-
function (x, max.dist, vario.col = "all", scaled = FALSE, var.lines = FALSE, 
    envelope.obj = NULL, pts.range.cex, bin.cloud = FALSE, ...) 
{
    if (missing(max.dist)) 
        max.dist <- max(x$u)
    Ldots <- list(...)
    if (is.null(Ldots$xlab)) 
        Ldots$xlab <- "distance"
    if (is.null(Ldots$ylab)) 
        Ldots$ylab <- "semi-variance"
    if (is.null(Ldots$ty)) {
        if (x$output.type == "bin") 
            Ldots$type <- "p"
        if (x$output.type == "smooth") 
            Ldots$type <- "l"
        if (x$output.type == "cloud") 
            Ldots$type <- "p"
    }
    if (is.null(Ldots$col)) 
        Ldots$col <- 1:6
    if (is.null(Ldots$lty)) 
        Ldots$lty <- 1:5
    if (is.null(Ldots$lwd)) 
        Ldots$lwd <- 1
    if (is.null(Ldots$pch)) 
        Ldots$pch <- NULL
    if (is.null(Ldots$cex)) 
        Ldots$cex <- NULL
    if (is.null(Ldots$add)) 
        Ldots$add <- FALSE
    if (bin.cloud == TRUE && all(is.na(x$bin.cloud))) 
        stop("plot.variogram: object must be a binned variogram with option bin.cloud=TRUE")
    if (bin.cloud == TRUE && any(!is.na(x$bin.cloud))) 
        boxplot(x$bin.cloud, varwidth = TRUE, xlab = Ldots$xlab, 
            ylab = Ldots$ylab,cex.lab=1.2,pch=3,cex=0.5)
    else {
        if (!missing(pts.range.cex)) {
            cex.min <- min(pts.range.cex)
            cex.max <- max(pts.range.cex)
            if (cex.min != cex.max) {
                pts.prop <- TRUE
                sqn <- sqrt(x$n[x$u <= max.dist])
                pts.cex <- cex.min + ((sqn - min(sqn)) * (cex.max - 
                  cex.min)/(max(sqn) - min(sqn)))
            }
            else pts.prop <- FALSE
        }
        else pts.prop <- FALSE
        u <- x$u[x$u <= max.dist]
        v <- x$v
        if (is.vector(v) | length(v) == length(x$u)) 
            v <- matrix(v, ncol = 1)
        v <- v[x$u <= max.dist, , drop = FALSE]
        if (vario.col == "all") 
            vario.col <- 1:dim(v)[2]
        else if (mode(vario.col) != "numeric" | any(vario.col > 
            ncol(v))) 
            stop("argument vario.col must be equals to \"all\" or a vector indicating the column numbers to be plotted")
        v <- v[, vario.col, drop = FALSE]
        if (scaled) 
            v <- t(t(v)/x$var.mark[vario.col])
        if (is.null(Ldots$ylim)) {
            ymax <- max(v)
            if (!is.null(envelope.obj)) 
                ymax <- max(c(envelope.obj$v.upper, ymax))
            Ldots$ylim <- c(0, ymax)
        }
        if (ncol(v) == 1) {
            v <- as.vector(v)
            uv <- data.frame(distance = u, semivariance = v)
            if (is.null(list(...)$ylim)) {
                if (pts.prop) 
                  plot(uv, xlim = c(0, max.dist), ylim = Ldots$ylim, 
                    cex = pts.cex, ...)
                else plot(uv, xlim = c(0, max.dist), ylim = Ldots$ylim, 
                  ...)
            }
            else {
                if (pts.prop) 
                  plot(uv, xlim = c(0, max.dist), ylim = Ldots$ylim, 
                    cex = pts.cex)
                else plot(uv, xlim = c(0, max.dist), ylim = Ldots$ylim)
            }
        }
        else matplot(x = u, y = v, xlim = c(0, max.dist), ylim = Ldots$ylim, 
            xlab = Ldots$xlab, ylab = Ldots$ylab, type = Ldots$type, 
            add = Ldots$add, pch = Ldots$pch, lty = Ldots$lty, 
            lwd = Ldots$lwd, col = Ldots$col)
        if (var.lines) {
            if (scaled) 
                abline(h = 1, lty = 3)
            else abline(h = x$var.mark, lty = 3)
        }
        if (!is.null(envelope.obj)) {
            lines(u, envelope.obj$v.lower, lty = 4)
            lines(u, envelope.obj$v.upper, lty = 4)
        }
    }
    return(invisible())
}
