ecdfPlot <-
function (x, discrete = FALSE, prob.method = ifelse(discrete, 
    "emp.probs", "plot.pos"), plot.pos.con = 0.375, plot.it = TRUE, 
    add = FALSE, ecdf.col = "black", ecdf.lwd = 3 * par("cex"), 
    ecdf.lty = 1, curve.fill = FALSE, curve.fill.col = "cyan", 
    ..., type = ifelse(discrete, "s", "l"), main = NULL, xlab = NULL, 
    ylab = NULL, xlim = NULL, ylim = NULL) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if (!is.vector(plot.pos.con, mode = "numeric") || length(plot.pos.con) != 
        1 || plot.pos.con < 0 || plot.pos.con > 1) 
        stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    n <- length(x)
    ord.stats <- sort(x)
    prob.method <- match.arg(prob.method, c("emp.probs", "plot.pos"))
    if (prob.method == "emp.probs") {
        rle.list <- rle(ord.stats)
        F.x <- cumsum(rle.list$lengths)/n
        ord.stats <- rle.list$values
    }
    else F.x <- ppoints(n, a = plot.pos.con)
    if (plot.it) {
        if (!add) {
            if (is.null(main)) 
                main <- paste("Empirical CDF of", data.name)
            if (is.null(xlab)) 
                xlab <- paste("Order Statistics for", data.name)
            if (is.null(ylab)) 
                ylab <- "Cumulative Probability"
            if (is.null(xlim)) 
                xlim <- range(ord.stats)
            if (is.null(ylim)) 
                ylim <- c(0, 1)
            plot(ord.stats, F.x, type = "n", ..., xlim = xlim, 
                ylim = ylim, xlab = xlab, ylab = ylab, main = main)
            arg.list <- list(x = ord.stats, y = F.x)
            arg.list <- c(arg.list, checkGraphicsPars(...)$gen.gp.list, 
                list(type = type, col = ecdf.col, lwd = ecdf.lwd, 
                  lty = ecdf.lty))
            do.call("lines", arg.list)
        }
        else lines(ord.stats, F.x, ..., type = type, col = ecdf.col, 
            lwd = ecdf.lwd, lty = ecdf.lty)
        if (prob.method == "emp.probs") 
            lines(ord.stats[c(1, 1)], c(0, F.x[1]), ..., type = type, 
                col = ecdf.col, lwd = ecdf.lwd, lty = ecdf.lty)
        if (curve.fill) {
            n <- length(F.x)
            if (type != "s") {
                polygon(c(ord.stats, rev(ord.stats)), c(F.x, 
                  rep(0, n)), border = FALSE, col = curve.fill.col)
            }
            else {
                dum.x <- as.vector(matrix(ord.stats, 2, n, byrow = TRUE))
                m <- 2 * n
                dum.x <- dum.x[-c(1, m)]
                dum.y <- as.vector(matrix(F.x, 2, n, byrow = TRUE))
                dum.y <- dum.y[-c(m - 1, m)]
                polygon(c(dum.x, rev(dum.x)), c(dum.y, rep(0, 
                  m - 2)), border = FALSE, col = curve.fill.col)
            }
        }
    }
    invisible(list(Order.Statistics = ord.stats, Cumulative.Probabilities = F.x))
}
