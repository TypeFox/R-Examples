ecdfPlotCensored <-
function (x, censored, censoring.side = "left", discrete = FALSE, 
    prob.method = "michael-schucany", plot.pos.con = 0.375, plot.it = TRUE, 
    add = FALSE, ecdf.col = 1, ecdf.lwd = 3 * par("cex"), ecdf.lty = 1, 
    include.cen = FALSE, cen.pch = ifelse(censoring.side == "left", 
        6, 2), cen.cex = par("cex"), cen.col = 4, ..., type = ifelse(discrete, 
        "s", "l"), main = NULL, xlab = NULL, ylab = NULL, xlim = NULL, 
    ylim = NULL) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    if (!is.vector(censored, mode = "numeric") & !is.vector(censored, 
        mode = "logical")) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    data.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(as.numeric(censored))))) > 
        0) {
        is.not.finite.warning(x)
        is.not.finite.warning(as.numeric(censored))
        x <- x[ok]
        censored <- censored[ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and 'censored' removed."))
    }
    if (is.numeric(censored)) {
        if (!all(censored == 0 | censored == 1)) 
            stop(paste("When 'censored' is a numeric vector, all values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
        censored <- as.logical(censored)
    }
    if (sum(censored) == 0) 
        stop(paste("No censored values indicated by 'censored';", 
            "use \n\t\t\tthe function 'ecdfPlot'"))
    if (length(unique(x[!censored])) < 1) 
        stop("'x' must contain at least one non-missing, uncensored value.")
    prob.method <- match.arg(prob.method, c("michael-schucany", 
        "hirsch-stedinger", "kaplan-meier", "nelson"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (censoring.side == "left" && prob.method == "nelson") 
        stop("Nelson method not available for censoring.side='left'")
    if (!is.vector(plot.pos.con, mode = "numeric") || length(plot.pos.con) != 
        1 || plot.pos.con < 0 || plot.pos.con > 1) 
        stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
    ppoints.list <- ppointsCensored(x = x, censored = censored, 
        censoring.side = censoring.side, prob.method = prob.method, 
        plot.pos.con = plot.pos.con)
    x <- ppoints.list$Order.Statistics
    p <- ppoints.list$Cumulative.Probabilities
    cen <- ppoints.list$Censored
    if (plot.it) {
        x.no.cen <- x[!cen]
        p.no.cen <- p[!cen]
        if (!add) {
            if (is.null(main)) {
                prob.method.name = switch(prob.method, `michael-schucany` = "Michael-Schucany", 
                  `hirsch-stedinger` = "Hirsch-Stedinger", `kaplan-meier` = "Kaplan-Meier", 
                  nelson = "Nelson")
                main <- paste("Empirical CDF of ", data.name, 
                  ", Based on\n", prob.method.name, " Plotting Positions (Censored Data)", 
                  sep = "")
            }
            if (is.null(xlab)) 
                xlab <- paste("Order Statistics for", data.name)
            if (is.null(ylab)) 
                ylab <- "Cumulative Probability"
            if (is.null(xlim)) 
                xlim <- range(x)
            if (is.null(ylim)) 
                ylim <- c(0, 1)
            plot(x.no.cen, p.no.cen, type = "n", ..., xlim = xlim, 
                ylim = ylim, xlab = xlab, ylab = ylab, main = main)
            arg.list <- list(x = x.no.cen, y = p.no.cen)
            arg.list <- c(arg.list, checkGraphicsPars(...)$gen.gp.list, 
                list(type = type, col = ecdf.col, lwd = ecdf.lwd, 
                  lty = ecdf.lty))
            do.call("lines", arg.list)
        }
        else lines(x.no.cen, p.no.cen, ..., type = type, col = ecdf.col, 
            lwd = ecdf.lwd, lty = ecdf.lty)
        if (type == "s") 
            lines(x.no.cen[c(1, 1)], c(0, p.no.cen[1]), ..., 
                type = type, col = ecdf.col, lwd = ecdf.lwd, 
                lty = ecdf.lty)
    }
    if (include.cen) 
        points(x[cen], p[cen], pch = cen.pch, cex = cen.cex, 
            col = cen.col)
    invisible(ppoints.list)
}
