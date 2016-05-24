cdfCompare <-
function (x, y = NULL, discrete = FALSE, prob.method = ifelse(discrete, 
    "emp.probs", "plot.pos"), plot.pos.con = NULL, distribution = "norm", 
    param.list = NULL, estimate.params = is.null(param.list), 
    est.arg.list = NULL, x.col = "blue", y.or.fitted.col = "black", 
    x.lwd = 3 * par("cex"), y.or.fitted.lwd = 3 * par("cex"), 
    x.lty = 1, y.or.fitted.lty = 2, digits = .Options$digits, 
    ..., type = ifelse(discrete, "s", "l"), main = NULL, xlab = NULL, 
    ylab = NULL, xlim = NULL, ylim = NULL) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    x.name <- deparse(substitute(x))
    if ((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
        is.not.finite.warning(x)
        x <- x[x.ok]
        warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    check.gp.list <- checkGraphicsPars(...)
    gp.arg.list <- check.gp.list$gp.arg.list
    gen.gp.list <- check.gp.list$gen.gp.list
    if (!is.null(y)) {
        if (!is.vector(y, mode = "numeric") || is.factor(y)) 
            stop("'y' must be a numeric vector")
        y.name <- deparse(substitute(y))
        if ((bad.obs <- sum(!(y.ok <- is.finite(y)))) > 0) {
            is.not.finite.warning(y)
            y <- y[y.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'y' removed."))
        }
        if (is.null(plot.pos.con)) {
            plot.pos.con <- 0.375
        }
        else {
            if (!is.vector(plot.pos.con, mode = "numeric") || 
                length(plot.pos.con) != 1 || plot.pos.con < 0 || 
                plot.pos.con > 1) 
                stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
        }
        if (is.null(main)) 
            main <- paste("Empirical CDF for ", x.name, " (solid line)\nwith Empirical CDF for ", 
                y.name, " (dashed line)", sep = "", collapse = "")
        if (is.null(xlab)) 
            xlab <- paste("Order Statistics for", x.name, "and", 
                y.name)
        if (is.null(xlim)) 
            xlim <- range(x, y)
        x.ecdf.list <- ecdfPlot(x, discrete = discrete, prob.method = prob.method, 
            plot.pos.con = plot.pos.con, ecdf.col = x.col, ecdf.lwd = x.lwd, 
            ecdf.lty = x.lty, ..., type = type, main = main, 
            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
        arg.list <- list(x = y, discrete = discrete, prob.method = prob.method, 
            plot.pos.con = plot.pos.con, add = TRUE, ecdf.col = y.or.fitted.col, 
            ecdf.lwd = y.or.fitted.lwd, ecdf.lty = y.or.fitted.lty)
        arg.list <- c(arg.list, gen.gp.list, list(type = type))
        y.ecdf.list <- do.call("ecdfPlot", arg.list)
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
            if (EnvStats::Distribution.df[dist.abb, "Estimation.Method(s)"] == 
                "") 
                stop(paste("No estimation method available for the", 
                  dist.name, "Distribution"))
            ename <- paste("e", dist.abb, sep = "")
            est.list <- do.call(ename, c(list(x = x), est.arg.list))
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
            main <- paste("Empirical CDF for ", x.name, " (solid line)\nwith Fitted ", 
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
        x.ecdf.list <- ecdfPlot(x, discrete = discrete, prob.method = prob.method, 
            plot.pos.con = plot.pos.con, ecdf.col = x.col, ecdf.lwd = x.lwd, 
            ecdf.lty = x.lty, ..., type = type, main = main, 
            xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
        arg.list <- list(distribution = dist.abb, param.list = param.list, 
            add = TRUE, cdf.col = y.or.fitted.col, cdf.lwd = y.or.fitted.lwd, 
            cdf.lty = y.or.fitted.lty)
        arg.list <- c(arg.list, gen.gp.list, list(type = type))
        fitted.cdf.list <- do.call("cdfPlot", arg.list)
        ret.list <- list(x.ecdf.list = x.ecdf.list, fitted.cdf.list = fitted.cdf.list)
    }
    invisible(ret.list)
}
