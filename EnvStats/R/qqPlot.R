qqPlot <-
function (x, y = NULL, distribution = "norm", param.list = list(mean = 0, 
    sd = 1), estimate.params = plot.type == "Tukey Mean-Difference Q-Q", 
    est.arg.list = NULL, plot.type = "Q-Q", plot.pos.con = NULL, 
    plot.it = TRUE, equal.axes = qq.line.type == "0-1" || estimate.params, 
    add.line = FALSE, qq.line.type = "least squares", duplicate.points.method = "standard", 
    points.col = 1, line.col = 1, line.lwd = par("cex"), line.lty = 1, 
    digits = .Options$digits, ..., main = NULL, xlab = NULL, 
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
    nx <- length(x)
    plot.type <- match.arg(plot.type, c("Q-Q", "Tukey Mean-Difference Q-Q"))
    duplicate.points.method <- match.arg(duplicate.points.method, 
        c("standard", "jitter", "number"))
    qq.line.type <- match.arg(qq.line.type, c("least squares", 
        "0-1", "robust"))
    gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    ret.list <- list()
    if (!is.null(y)) {
        if (!is.vector(y, mode = "numeric")) 
            stop("'y' must be a numeric vector")
        y.name <- deparse(substitute(y))
        if ((bad.obs <- sum(!(y.ok <- is.finite(y)))) > 0) {
            is.not.finite.warning(y)
            y <- y[y.ok]
            warning(paste(bad.obs, "observations with NA/NaN/Inf in 'y' removed."))
        }
        if (is.null(plot.pos.con)) {
            plot.pos.con <- 0.5
        }
        else {
            if (!is.vector(plot.pos.con, mode = "numeric") || 
                length(plot.pos.con) != 1 || plot.pos.con < 0 || 
                plot.pos.con > 1) 
                stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
        }
        q.x <- sort(x)
        q.y <- sort(y)
        ny <- length(y)
        if (nx != ny) {
            ppoints.x <- ppoints(nx, a = plot.pos.con)
            ppoints.y <- ppoints(ny, a = plot.pos.con)
            if (nx > ny) {
                q.x <- approx(ppoints.x, q.x, xout = ppoints.y, 
                  rule = 2)$y
            }
            else {
                q.y <- approx(ppoints.y, q.y, xout = ppoints.x, 
                  rule = 2)$y
            }
        }
        if (plot.it) {
            if (plot.type == "Q-Q") {
                if (is.null(xlab)) 
                  xlab <- paste("Quantiles of", x.name)
                if (is.null(ylab)) 
                  ylab <- paste("Quantiles of", y.name)
                if (is.null(main)) 
                  main <- paste("Q-Q Plot of\n", y.name, "vs.", 
                    x.name)
            }
            else {
                if (is.null(xlab)) 
                  xlab <- "Mean of Quantiles"
                if (is.null(ylab)) 
                  ylab <- paste(y.name, "Quantiles -", x.name, 
                    "Quantiles")
                if (is.null(main)) 
                  main <- paste("Tukey Mean-Difference Q-Q Plot for\n", 
                    x.name, "and", y.name)
            }
        }
    }
    else {
        check.da.list <- check.distribution.args(distribution, 
            check.params = FALSE)
        distribution <- check.da.list$dist.abb
        special <- any(distribution == c("lnorm", "lnorm3", "zmlnorm", 
            "zmnorm"))
        if (!estimate.params) {
            if (distribution != "norm" && !special && missing(param.list)) 
                stop(paste("When 'estimate.params=F' you must supply", 
                  "the argument 'param.list'"))
            check.da.list <- check.distribution.args(distribution, 
                param.list)
        }
        dist.name <- check.da.list$dist.name
        n.dist.params <- check.da.list$n.dist.params
        dist.params.names <- check.da.list$dist.params.names
        zm <- any(distribution == c("zmnorm", "zmlnorm", "zmlnormAlt"))
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
        if (estimate.params) {
            if (EnvStats::Distribution.df[distribution, "Estimation.Method(s)"] == 
                "") 
                stop(paste("No estimation method available for the", 
                  dist.name, "Distribution"))
            ename <- paste("e", distribution, sep = "")
            est.param.vec <- do.call(ename, c(list(x = x), est.arg.list))$parameters
            if (dist.params.names[n.dist.params] == "ncp") {
                warning(paste("No estimation method available for", 
                  "Non-Central Distributions.\n "))
                est.param.vec <- c(est.param.vec, 0)
            }
            if (zm) 
                est.param.vec <- est.param.vec[1:n.dist.params]
            param.list <- as.list(est.param.vec)
            names(param.list) <- dist.params.names
        }
        else param.list <- check.da.list$param.list
        distribution.x <- ifelse(special, "norm", distribution)
        dist.name.x <- ifelse(special, "Normal", dist.name)
        if (distribution == "zmlnormAlt") {
            distribution.x <- "lnormAlt"
            dist.name.x <- "Lognormal"
        }
        q.y <- switch(distribution, lnorm = sort(log(x)), lnorm3 = sort(log(x - 
            param.list$threshold)), zmlnorm = sort(log(x[x > 
            0])), zmnorm = , zmlnormAlt = sort(x[x != 0]), sort(x))
        if (zm) 
            nx <- length(q.y)
        param.list.x <- param.list
        if (special) {
            switch(distribution, lnorm = {
            }, lnorm3 = , zmlnorm = param.list.x <- param.list[c("meanlog", 
                "sdlog")], zmnorm = param.list.x <- param.list[c("mean", 
                "sd")])
            names(param.list.x) <- c("mean", "sd")
        }
        if (distribution == "zmlnormAlt") {
            param.list.x <- param.list[c("mean", "cv")]
            names(param.list.x) <- c("mean", "cv")
        }
        qname <- paste("q", distribution.x, sep = "")
        q.x <- do.call(qname, c(list(ppoints(nx, a = plot.pos.con)), 
            param.list.x))
        dist.params.names.x <- names(param.list.x)
        if (plot.it) {
            qlab <- switch(distribution, lnorm = paste("Log[", 
                x.name, "]", sep = ""), lnorm3 = paste("Log[", 
                x.name, "-", format(param.list$threshold, digits = digits), 
                "]"), zmnorm = paste("Non-Zero Values of", x.name), 
                zmlnorm = paste("Log [", x.name, "> 0 ]"), zmlnormAlt = paste(x.name, 
                  "> 0"), x.name)
            if (any(distribution == c("beta", "chisq", "f")) && 
                param.list$ncp == 0) {
                x.string <- paste(dist.name, "(", paste(paste(dist.params.names[-n.dist.params], 
                  signif(unlist(param.list[-n.dist.params]), 
                    digits), sep = " = "), collapse = ", "), 
                  ")", sep = "")
            }
            else {
                x.string <- paste(dist.name.x, "(", paste(paste(dist.params.names.x, 
                  signif(unlist(param.list.x), digits), sep = " = "), 
                  collapse = ", "), ")", sep = "")
            }
            if (plot.type == "Q-Q") {
                if (is.null(xlab)) {
                  xlab <- paste("Quantiles of", x.string)
                }
                if (is.null(ylab)) 
                  ylab <- paste("Quantiles of", qlab)
                if (is.null(main)) 
                  main <- paste(dist.name.x, "Q-Q Plot for", 
                    qlab)
            }
            else {
                if (is.null(xlab)) 
                  xlab <- "Mean of Observed and Fitted Quantiles"
                if (is.null(ylab)) 
                  ylab <- "Observed - Fitted Quantiles"
                if (is.null(main)) 
                  main <- paste("Tukey Mean-Difference Q-Q Plot for ", 
                    qlab, "\nFitted to ", x.string, " Distribution", 
                    sep = "")
            }
        }
    }
    if (plot.type == "Q-Q") {
        if (plot.it) {
            if (is.null(xlim) && is.null(ylim) && equal.axes) {
                xlim <- range(q.x, q.y)
                ylim <- xlim
            }
            else {
                if (is.null(xlim)) 
                  xlim <- range(q.x)
                if (is.null(ylim)) 
                  ylim <- range(q.y)
            }
            plot(q.x, q.y, type = "n", ..., xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, main = main)
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
        ret.list <- list(x = q.x, y = q.y)
    }
    else {
        q.mean <- (q.x + q.y)/2
        q.diff <- q.y - q.x
        if (plot.it) {
            if (is.null(ylim)) {
                rqmo2 <- diff(range(q.mean))/2
                mqd <- median(q.diff)
                ylim.min <- min(min(q.diff), mqd - rqmo2)
                ylim.max <- max(max(q.diff), mqd + rqmo2)
                ylim <- c(ylim.min, ylim.max)
            }
            if (is.null(xlim)) 
                xlim <- range(q.mean)
            plot(q.mean, q.diff, type = "n", ..., xlim = xlim, 
                ylim = ylim, xlab = xlab, ylab = ylab, main = main)
            arg.list <- c(list(x = q.mean, y = q.diff, method = duplicate.points.method), 
                gen.gp.list, list(col = points.col))
            do.call("points.w.dups", arg.list)
            if (add.line) {
                arg.list <- c(list(h = 0), gen.gp.list, list(col = line.col, 
                  lwd = line.lwd, lty = line.lty))
                do.call("abline", arg.list)
            }
        }
        ret.list <- list(x = q.mean, y = q.diff)
    }
    invisible(ret.list)
}
