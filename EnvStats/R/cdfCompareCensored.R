cdfCompareCensored <-
function (x, censored, censoring.side = "left", y = NULL, y.censored = NULL, 
    y.censoring.side = censoring.side, discrete = FALSE, prob.method = "michael-schucany", 
    plot.pos.con = NULL, distribution = "norm", param.list = NULL, 
    estimate.params = is.null(param.list), est.arg.list = NULL, 
    x.col = "blue", y.or.fitted.col = "black", x.lwd = 3 * par("cex"), 
    y.or.fitted.lwd = 3 * par("cex"), x.lty = 1, y.or.fitted.lty = 2, 
    include.x.cen = FALSE, x.cen.pch = ifelse(censoring.side == 
        "left", 6, 2), x.cen.cex = par("cex"), x.cen.col = "red", 
    include.y.cen = FALSE, y.cen.pch = ifelse(y.censoring.side == 
        "left", 6, 2), y.cen.cex = par("cex"), y.cen.col = "black", 
    digits = .Options$digits, ..., type = ifelse(discrete, "s", 
        "l"), main = NULL, xlab = NULL, ylab = NULL, xlim = NULL, 
    ylim = NULL) 
{
    if (!is.vector(x, mode = "numeric")) 
        stop("'x' must be a numeric vector")
    if (!is.vector(censored, mode = "numeric") & !is.vector(censored, 
        mode = "logical")) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    x.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(ok <- is.finite(x) & is.finite(as.numeric(censored))))) > 
        0) {
        x <- x[ok]
        censored <- censored[ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' and/or 'censored' removed."))
    }
    if (is.numeric(censored)) {
        if (!all(censored == 0 | censored == 1)) 
            stop(paste("When 'censored' is a numeric vector, all values of", 
                "'censored' must be 0 (not censored) or 1 (censored)."))
        censored <- as.logical(censored)
    }
    n.cen.x <- sum(censored)
    K <- sort(unique(x[censored]))
    if (n.cen.x == 0 && is.null(y)) 
        stop(paste("No censored values indicated by 'censored';", 
            "use \n\t\t\tthe function 'cdfCompare'"))
    if (length(unique(x[!censored])) < 2) 
        stop("'x' must contain at least two non-missing, uncensored value.")
    check.gp.list <- checkGraphicsPars(...)
    gp.arg.list <- check.gp.list$gp.arg.list
    gen.gp.list <- check.gp.list$gen.gp.list
    prob.method <- match.arg(prob.method, c("michael-schucany", 
        "hirsch-stedinger", "kaplan-meier", "nelson"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (n.cen.x > 0 && censoring.side == "left" && prob.method == 
        "nelson") 
        stop("Nelson method not available for censoring.side='left'")
    if (!is.null(y)) {
        if (!is.vector(y, mode = "numeric")) 
            stop("'y' must be a numeric vector")
        if (!is.vector(y.censored, mode = "numeric") & !is.vector(y.censored, 
            mode = "logical")) 
            stop("'y.censored' must be a logical or numeric vector")
        if (length(y.censored) != length(y)) 
            stop("'y.censored' must be the same length as 'y'")
        y.name <- deparse(substitute(y))
        if ((bad.obs <- sum(!(ok <- is.finite(y) & is.finite(as.numeric(y.censored))))) > 
            0) {
            y <- y[ok]
            y.censored <- y.censored[ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'y' and/or 'y.censored' removed."))
        }
        if (is.numeric(y.censored)) {
            if (!all(y.censored == 0 | y.censored == 1)) 
                stop(paste("When 'y.censored' is a numeric vector, all values of", 
                  "'y.censored' must be 0 (not censored) or 1 (censored)."))
            y.censored <- as.logical(y.censored)
        }
        n.cen.y <- sum(y.censored)
        if (n.cen.x == 0 && n.cen.y == 0) 
            stop(paste("No censored values indicated by 'censored'", 
                "and 'y.censored';", "use \n\t\t\tthe function 'cdfCompare'"))
        if (length(unique(y[!y.censored])) < 2) 
            stop("'y' must contain at least two non-missing, uncensored value.")
        if (is.null(plot.pos.con)) {
            plot.pos.con <- 0.375
        }
        else {
            if (!is.vector(plot.pos.con, mode = "numeric") || 
                length(plot.pos.con) != 1 || plot.pos.con < 0 || 
                plot.pos.con > 1) 
                stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
        }
        y.censoring.side <- match.arg(y.censoring.side, c("left", 
            "right"))
        if (n.cen.y > 0 && y.censoring.side == "left" && prob.method == 
            "nelson") 
            stop("Nelson method not available for y.censoring.side='left'")
        if (is.null(main)) {
            if (n.cen.x > 0) 
                x.string <- paste(x.name, "(Censored; solid line)\n")
            else x.string <- paste(x.name, "(solid line)\n")
            if (n.cen.y > 0) 
                y.string <- paste(y.name, "(Censored; dashed line)")
            else y.string <- paste(y.name, "(dashed line)\n")
            main <- paste("Empirical CDF for ", x.string, "with Empirical CDF for ", 
                y.string, sep = "", collapse = "")
        }
        if (is.null(xlab)) 
            xlab <- paste("Order Statistics for", x.name, "and", 
                y.name)
        if (is.null(xlim)) 
            xlim <- range(x, y)
        if (n.cen.x == 0) 
            x.ecdf.list <- ecdfPlot(x, discrete = discrete, prob.method = "plot.pos", 
                plot.pos.con = plot.pos.con, ecdf.col = x.col, 
                ecdf.lwd = x.lwd, ecdf.lty = x.lty, ..., type = type, 
                main = main, xlab = xlab, ylab = ylab, xlim = xlim, 
                ylim = ylim)
        else x.ecdf.list <- ecdfPlotCensored(x = x, censored = censored, 
            censoring.side = censoring.side, discrete = discrete, 
            prob.method = prob.method, plot.pos.con = plot.pos.con, 
            ecdf.col = x.col, ecdf.lwd = x.lwd, ecdf.lty = x.lty, 
            include.cen = include.x.cen, cen.pch = x.cen.pch, 
            cen.cex = x.cen.cex, cen.col = x.cen.col, ..., type = type, 
            main = main, xlab = xlab, ylab = ylab, xlim = xlim, 
            ylim = ylim)
        if (n.cen.y == 0) {
            arg.list <- list(x = y, discrete = discrete, prob.method = "plot.pos", 
                plot.pos.con = plot.pos.con, add = TRUE, ecdf.col = y.or.fitted.col, 
                ecdf.lwd = y.or.fitted.lwd, ecdf.lty = y.or.fitted.lty)
            arg.list <- c(arg.list, gen.gp.list, list(type = type))
            y.ecdf.list <- do.call("ecdfPlot", arg.list)
        }
        else {
            arg.list <- list(x = y, censored = y.censored, censoring.side = y.censoring.side, 
                discrete = discrete, prob.method = prob.method, 
                plot.pos.con = plot.pos.con, add = TRUE, ecdf.col = y.or.fitted.col, 
                ecdf.lwd = y.or.fitted.lwd, ecdf.lty = y.or.fitted.lty, 
                include.cen = include.y.cen, cen.pch = y.cen.pch, 
                cen.cex = y.cen.cex, cen.col = y.cen.col)
            arg.list <- c(arg.list, gen.gp.list, list(type = type))
            y.ecdf.list <- do.call("ecdfPlotCensored", arg.list)
        }
        ret.list <- list(x.ecdf.list = x.ecdf.list, y.ecdf.list = y.ecdf.list)
    }
    else {
        if (estimate.params) 
            check.da.list <- check.distribution.args(distribution, 
                check.params = FALSE)
        else {
            if (is.null(param.list)) 
                stop(paste("When 'estimate.params=F' you must supply", 
                  "the argument 'param.list'"))
            check.da.list <- check.distribution.args(distribution, 
                param.list)
        }
        dist.abb <- check.da.list$dist.abb
        dist.name <- check.da.list$dist.name
        dist.type <- check.da.list$dist.type
        n.dist.params <- check.da.list$n.dist.params
        dist.params.names <- check.da.list$dist.params.names
        if (estimate.params) {
            e.fcn <- paste("e", dist.abb, "Censored", sep = "")
            if (length(find(e.fcn)) == 0) 
                stop(paste("No estimation method available for the", 
                  dist.name, "Distribution with Censored Data"))
            est.list <- do.call(e.fcn, c(list(x = x, censored = censored, 
                censoring.side = censoring.side), est.arg.list))
            est.param.vec <- est.list$parameters
            estimation.method <- est.list$method
            if (dist.params.names[n.dist.params] == "ncp") {
                warning(paste("No estimation method available for", 
                  "Non-Central Distributions.\n "))
                est.param.vec <- c(est.param.vec, 0)
            }
            param.list <- as.list(est.param.vec)
            names(param.list) <- dist.params.names
        }
        else param.list <- check.da.list$param.list
        if (is.null(plot.pos.con)) {
            plot.pos.con <- switch(distribution, norm = , lnorm = , 
                lnormAlt = , lnorm3 = , zmnorm = , zmlnorm = , 
                zmlnormAlt = 0.375, evd = 0.44, 0.4)
        }
        else {
            if (!is.vector(plot.pos.con, mode = "numeric") || 
                length(plot.pos.con) != 1 || plot.pos.con < 0 || 
                plot.pos.con > 1) 
                stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
        }
        discrete <- dist.type == "Discrete" || dist.type == "Finite Discrete"
        if (is.null(main)) {
            main <- paste("Empirical CDF for ", x.name, " (Censored; solid line)\nwith Fitted ", 
                dist.name, " CDF (dashed line)", sep = "", collapse = "")
        }
        if (is.null(xlab)) {
            if (any(dist.abb == c("beta", "chisq", "f", "t"))) {
                if (param.list$ncp > 0) 
                  string <- paste("Non-central ", dist.name, 
                    "(", paste(paste(dist.params.names, signif(unlist(param.list), 
                      digits), sep = "="), collapse = ", "), 
                    ")", sep = "")
                else {
                  string <- paste(dist.name, "(", paste(paste(dist.params.names[-n.dist.params], 
                    signif(unlist(param.list[-n.dist.params]), 
                      digits), sep = "="), collapse = ", "), 
                    ")", sep = "")
                }
            }
            else {
                string <- paste(dist.name, "(", paste(paste(dist.params.names, 
                  signif(unlist(param.list), digits), sep = "="), 
                  collapse = ", "), ")", sep = "")
            }
            xlab <- paste("Order Statistics for ", x.name, " and\n", 
                string, " Distribution", sep = "")
        }
        if (is.null(xlim)) {
            qname <- paste("q", dist.abb, sep = "")
            xlim <- do.call(qname, c(list(p = c(0.001, 0.999)), 
                param.list))
        }
        x.ecdf.list <- ecdfPlotCensored(x = x, censored = censored, 
            censoring.side = censoring.side, discrete = discrete, 
            prob.method = prob.method, plot.pos.con = plot.pos.con, 
            ecdf.col = x.col, ecdf.lwd = x.lwd, ecdf.lty = x.lty, 
            include.cen = include.x.cen, cen.pch = x.cen.pch, 
            cen.cex = x.cen.cex, cen.col = x.cen.col, ..., type = type, 
            main = main, xlab = xlab, ylab = ylab, xlim = xlim, 
            ylim = ylim)
        arg.list <- list(distribution = dist.abb, param.list = param.list, 
            add = TRUE, cdf.col = y.or.fitted.col, cdf.lwd = y.or.fitted.lwd, 
            cdf.lty = y.or.fitted.lty)
        arg.list <- c(arg.list, gen.gp.list, list(type = type))
        fitted.cdf.list <- do.call("cdfPlot", arg.list)
        ret.list <- list(x.ecdf.list = x.ecdf.list, fitted.cdf.list = fitted.cdf.list)
    }
    invisible(ret.list)
}
