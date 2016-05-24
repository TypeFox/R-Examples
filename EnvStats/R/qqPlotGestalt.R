qqPlotGestalt <-
function (distribution = "norm", param.list = list(mean = 0, 
    sd = 1), estimate.params = FALSE, est.arg.list = NULL, sample.size = 10, 
    num.pages = 2, num.plots.per.page = 4, nrow = ceiling(num.plots.per.page/2), 
    plot.type = "Q-Q", plot.pos.con = switch(dist.abb, norm = , 
        lnorm = , lnormAlt = , lnorm3 = 0.375, evd = 0.44, 0.4), 
    equal.axes = (qq.line.type == "0-1" || estimate.params), 
    margin.title = NULL, add.line = FALSE, qq.line.type = "least squares", 
    duplicate.points.method = "standard", points.col = 1, line.col = 1, 
    line.lwd = par("cex"), line.lty = 1, digits = .Options$digits, 
    same.window = TRUE, ask = same.window & num.pages > 1, mfrow = c(nrow, 
        num.plots.per.page/nrow), mar = c(4, 4, 1, 1) + 0.1, 
    oma = c(0, 0, 7, 0), mgp = c(2, 0.5, 0), ..., main = NULL, 
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL) 
{
    if (!((num.plots.per.page%%2) == 0 || num.plots.per.page == 
        1)) 
        stop("'num.plots.per.page' must be 1 or an even number")
    if ((num.plots.per.page%%nrow) != 0) 
        stop("'num.plots.per.page' must be a multiple of 'nrow'")
    plot.type <- match.arg(plot.type, c("Q-Q", "Tukey Mean-Difference Q-Q"))
    duplicate.points.method <- match.arg(duplicate.points.method, 
        c("standard", "jitter", "number"))
    qq.line.type <- match.arg(qq.line.type, c("least squares", 
        "0-1", "robust"))
    par.list <- list(mfrow = mfrow, mar = mar, oma = oma, mgp = mgp)
    dev.new()
    cex.orig <- par("cex")
    par(par.list)
    par(cex = 0.75 * cex.orig, mex = 0.75 * cex.orig)
    devAskNewPage(ask = ask)
    check.gp.list <- checkGraphicsPars(...)
    gp.arg.list <- check.gp.list$gp.arg.list
    gp.names <- check.gp.list$gp.names
    gen.gp.list <- check.gp.list$gen.gp.list
    check.da.list <- check.distribution.args(distribution, param.list)
    dist.abb <- check.da.list$dist.abb
    dist.name <- check.da.list$dist.name
    n.dist.params <- check.da.list$n.dist.params
    dist.params.names <- check.da.list$dist.params.names
    param.list <- check.da.list$param.list
    param.list.x <- param.list
    if (!is.vector(plot.pos.con, mode = "numeric") || length(plot.pos.con) != 
        1 || plot.pos.con < 0 || plot.pos.con > 1) 
        stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
    if (estimate.params) {
        if (EnvStats::Distribution.df[dist.abb, "Estimation.Method(s)"] == 
            "") 
            stop(paste("No estimation method available for the", 
                dist.name, "Distribution"))
        if (dist.params.names[n.dist.params] == "ncp" && param.list$ncp != 
            0) 
            stop("No estimation method available for Non-Central Distributions.")
    }
    r.fcn <- paste("r", dist.abb, sep = "")
    q.fcn <- paste("q", dist.abb, sep = "")
    if (is.null(margin.title)) {
        margin.title.supplied <- FALSE
        plot.string <- ifelse(num.plots.per.page == 1, "Plot", 
            "Plots")
        if (any(dist.abb == c("beta", "chisq", "f")) && param.list$ncp == 
            0) 
            margin.title <- paste(plot.type, " ", plot.string, 
                " for\n", dist.name, "(", paste(paste(dist.params.names[-n.dist.params], 
                  signif(unlist(param.list[-n.dist.params]), 
                    digits), sep = " = "), collapse = ", "), 
                ") Distribution", sep = "")
        else margin.title <- paste(plot.type, " ", plot.string, 
            " for\n", dist.name, "(", paste(paste(dist.params.names, 
                signif(unlist(param.list), digits), sep = " = "), 
                collapse = ", "), ") Distribution", sep = "")
    }
    else margin.title.supplied <- TRUE
    if (estimate.params) {
        est.fcn <- paste("e", dist.abb, sep = "")
        estimation.method <- do.call(est.fcn, c(list(x = do.call(r.fcn, 
            c(list(n = sample.size), param.list))), est.arg.list))$method
        ss.title <- paste("(Sample Size = ", sample.size, "; Estimation Method = ", 
            estimation.method, ")", sep = "")
    }
    else ss.title <- paste("(Sample Size = ", sample.size, "; No Parameter Estimation", 
        ")", sep = "")
    if (is.null(main)) 
        main <- ""
    if (plot.type == "Q-Q") {
        if (is.null(xlab)) 
            xlab <- "Quantiles of Assumed Distribution"
        if (is.null(ylab)) 
            ylab <- "Random Quantiles"
    }
    else {
        if (is.null(xlab)) 
            xlab <- "Mean of Observed and Fitted Quantiles"
        if (is.null(ylab)) 
            ylab <- "Observed-Fitted Quantiles"
    }
    user.xlim <- xlim
    user.ylim <- ylim
    for (i in 1:(num.pages * num.plots.per.page)) {
        if (i > 1 && ((i%%num.plots.per.page) == 1) && !same.window) {
            dev.new()
            par(par.list)
            par(cex = 0.75 * cex.orig, mex = 0.75 * cex.orig)
        }
        q.y <- sort(do.call(r.fcn, c(list(n = sample.size), param.list)))
        if (estimate.params) {
            est.param.vec <- do.call(est.fcn, c(list(x = q.y), 
                est.arg.list))$parameters
            if (!is.null(param.list$ncp)) 
                param.list.x[-n.dist.params] <- est.param.vec
            else param.list.x[] <- est.param.vec
        }
        q.x <- do.call(q.fcn, c(list(ppoints(q.y, a = plot.pos.con)), 
            param.list.x))
        if (plot.type == "Q-Q") {
            if (is.null(user.xlim) && is.null(user.ylim) && equal.axes) {
                xlim <- range(q.x, q.y)
                ylim <- xlim
            }
            else {
                if (is.null(user.xlim)) 
                  xlim <- range(q.x)
                if (is.null(user.ylim)) 
                  ylim <- range(q.y)
            }
            plot(q.x, q.y, type = "n", ..., main = main, xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim)
            arg.list <- c(list(x = q.x, y = q.y, method = duplicate.points.method), 
                gen.gp.list, list(col = points.col))
            do.call("points.w.dups", arg.list)
            if (add.line) 
                switch(qq.line.type, `least squares` = {
                  arg.list <- c(list(a = lm(q.y ~ q.x)), gen.gp.list, 
                    list(col = line.col, lwd = line.lwd, lty = line.lty))
                  do.call("abline", arg.list)
                }, `0-1` = {
                  arg.list <- c(list(a = 0, b = 1), gen.gp.list, 
                    list(col = line.col, lwd = line.lwd, lty = line.lty))
                  do.call("abline", arg.list)
                }, robust = {
                  arg.list <- c(list(x = q.x, y = q.y), gen.gp.list, 
                    list(col = line.col, lwd = line.lwd, lty = line.lty))
                  do.call("qqLine", arg.list)
                })
        }
        else {
            q.mean <- (q.x + q.y)/2
            q.diff <- q.y - q.x
            if (is.null(user.ylim)) {
                rqmo2 <- diff(range(q.mean))/2
                mqd <- median(q.diff)
                ylim.min <- min(min(q.diff), mqd - rqmo2)
                ylim.max <- max(max(q.diff), mqd + rqmo2)
                ylim <- c(ylim.min, ylim.max)
            }
            if (is.null(user.xlim)) 
                xlim <- range(q.mean)
            plot(q.mean, q.diff, type = "n", ..., main = main, 
                xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
            arg.list <- c(list(x = q.mean, y = q.diff, method = duplicate.points.method), 
                gen.gp.list, list(col = points.col))
            do.call("points.w.dups", arg.list)
            if (add.line) {
                arg.list <- c(list(h = 0), gen.gp.list, list(col = line.col, 
                  lwd = line.lwd, lty = line.lty))
                do.call("abline", arg.list)
            }
        }
        if ((i%%num.plots.per.page) == 0 && !margin.title.supplied) {
            mtext(margin.title, side = 3, line = 3, outer = TRUE, 
                cex = 1.25 * cex.orig)
            mtext(ss.title, side = 3, line = 0, outer = TRUE, 
                cex = cex.orig)
        }
    }
}
