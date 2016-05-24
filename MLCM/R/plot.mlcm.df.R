plot.mlcm.df <- function (x, clr = NULL, ...) 
{
	if(!is.factor(x$Resp)) x$Resp <- factor(x$Resp, levels = unique(x$Resp))
	# check to order so that all proportions map to upper left diagonal
	sw <- (x[, 2] > x[, 3]) | ((x[, 2] == x[, 3]) & (x[, 4] > x[, 5]))
	xl <- levels(x$Resp)
	x[sw, ] <- x[sw, c(1, 3:2, 5:4)]
	x$Resp <- unclass(x$Resp) - 1
	x[sw, 1] <- 1 - x[sw, 1]
	x$Resp <- factor(xl[x$Resp + 1], levels = xl)
    mx <- max(x[, -1])
    x$IntAct <- with(x, factor(x[, 2]):factor(x[, 4]):factor(x[, 
        3]):factor(x[, 5]))
    x.tab <- with(x, table(Resp, IntAct))
    x.prop <- 1 - apply(x.tab, 2, function(x) x/sum(x)) #x.tab/max(x.tab)
    x.mat <- matrix(x.prop[2, ], ncol = mx^2, byrow = TRUE)
    x.mat[col(x.mat) < row(x.mat)] <- NA
    clr <- if (is.null(clr)) {
        nr <- max(colSums(x.tab))
        stp <- trunc(83/nr)
        grey(seq(0, nr * stp, stp)/100)
    }
    else clr
    sl <- seq_len(mx)
    xy <- seq(0.5, mx + 0.5, len = mx^2 + 1)
    pmar <- par("mar")
    opar <- par(mar = pmar + c(1, 1, 0, 0), mgp = c(3.75, 1.75, 
        0))
    image(xy, xy, x.mat, col = clr, axes = FALSE, ...)
    box()
    abline(v = sl + 0.5, h = sl + 0.5)
    axis(1, at = sl, tick = FALSE, cex.axis = 2)
    axis(2, at = sl, tick = FALSE, cex.axis = 2)
    axis(3, at = 0:1 + 0.5, labels = range(sl), cex.axis = 0.9, 
        mgp = c(3, 0.75, 0))
    axis(4, at = 0:1 + 0.5, labels = range(sl), cex.axis = 0.9, 
        las = 2, mgp = c(3, 0.75, 0))
    par(opar)
    invisible()
}
