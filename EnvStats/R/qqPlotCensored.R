qqPlotCensored <-
function (x, censored, censoring.side = "left", prob.method = "michael-schucany", 
    plot.pos.con = NULL, distribution = "norm", param.list = list(mean = 0, 
        sd = 1), estimate.params = plot.type == "Tukey Mean-Difference Q-Q", 
    est.arg.list = NULL, plot.type = "Q-Q", plot.it = TRUE, equal.axes = qq.line.type == 
        "0-1" || estimate.params, add.line = FALSE, qq.line.type = "least squares", 
    duplicate.points.method = "standard", points.col = 1, line.col = 1, 
    line.lwd = par("cex"), line.lty = 1, digits = .Options$digits, 
    include.cen = FALSE, cen.pch = ifelse(censoring.side == "left", 
        6, 2), cen.cex = par("cex"), cen.col = 4, ..., main = NULL, 
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL) 
{
    if (!is.vector(x, mode = "numeric") || is.factor(x)) 
        stop("'x' must be a numeric vector")
    if (!((is.vector(censored, mode = "numeric") && !is.factor(censored)) || 
        is.vector(censored, mode = "logical"))) 
        stop("'censored' must be a logical or numeric vector")
    if (length(censored) != length(x)) 
        stop("'censored' must be the same length as 'x'")
    x.name <- deparse(substitute(x))
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
    cen.levels <- sort(unique(x[censored]))
    if (length(cen.levels) == 0) 
        stop(paste("No censored values indicated by 'censored';", 
            "use \n\t\t\tthe function 'qqPlot'"))
    x.no.cen <- x[!censored]
    if (length(unique(x.no.cen)) < 2) 
        stop("'x' must contain at least two non-missing, distinct uncensored value.")
    prob.method <- match.arg(prob.method, c("michael-schucany", 
        "hirsch-stedinger", "kaplan-meier", "modified kaplan-meier", 
        "nelson"))
    censoring.side <- match.arg(censoring.side, c("left", "right"))
    if (censoring.side == "left" && prob.method == "nelson") 
        stop("Nelson method not available for censoring.side='left'")
    if (censoring.side == "right" && prob.method == "modified kaplan-meier") 
        stop("Modified Kaplan-Meier method only relevant when censoring.side='left'")
    plot.type <- match.arg(plot.type, c("Q-Q", "Tukey Mean-Difference Q-Q"))
    duplicate.points.method <- match.arg(duplicate.points.method, 
        c("standard", "jitter", "number"))
    qq.line.type <- match.arg(qq.line.type, c("least squares", 
        "0-1", "robust"))
    gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    check.da.list <- check.distribution.args(distribution, check.params = FALSE)
    dist.name <- check.da.list$dist.name
    dist.abb <- check.da.list$dist.abb
    qname <- paste("q", dist.abb, sep = "")
    zm <- any(dist.abb == c("zmnorm", "zmlnorm", "zmlnormAlt"))
    log.dist <- any(dist.abb == c("lnorm", "lnorm3", "zmlnorm"))
    special <- log.dist || dist.abb == "zmnorm"
    if (!estimate.params) {
        if (distribution != "norm" && !special && missing(param.list)) 
            stop(paste("When 'estimate.params=F' you must supply", 
                "the argument 'param.list'"))
        check.da.list <- check.distribution.args(distribution, 
            param.list)
    }
    n.dist.params <- check.da.list$n.dist.params
    dist.params.names <- check.da.list$dist.params.names
    if (is.null(plot.pos.con)) {
        plot.pos.con <- switch(dist.abb, norm = , lnorm = , lnormAlt = , 
            lnorm3 = , zmnorm = , zmlnorm = , zmlnormAlt = 0.375, 
            evd = 0.44, 0.4)
    }
    else {
        if (!is.vector(plot.pos.con, mode = "numeric") || length(plot.pos.con) != 
            1 || plot.pos.con < 0 || plot.pos.con > 1) 
            stop("'plot.pos.con' must be a numeric scalar between 0 and 1")
    }
    if (estimate.params || (plot.it && plot.type == "Q-Q" && 
        is.null(ylim))) {
        e.fcn <- paste("e", dist.abb, "Censored", sep = "")
        if (length(find(e.fcn)) == 0) {
            if (estimate.params) 
                stop(paste("No estimation method available for the", 
                  dist.name, "Distribution with Censored Data"))
            else est.param.list <- NULL
        }
        else {
            est.param.vec <- do.call(e.fcn, c(list(x = x, censored = censored, 
                censoring.side = censoring.side), est.arg.list))$parameters
            if (dist.params.names[n.dist.params] == "ncp") {
                warning(paste("No estimation method available for", 
                  "Non-Central Distributions.\n "))
                est.param.vec <- c(est.param.vec, 0)
            }
            if (zm) 
                est.param.vec <- est.param.vec[1:n.dist.params]
            est.param.list <- as.list(est.param.vec)
            names(est.param.list) <- dist.params.names
            if (estimate.params) 
                param.list <- est.param.list
            else param.list <- check.da.list$param.list
        }
    }
    else param.list <- check.da.list$param.list
    dist.abb.x <- ifelse(special, "norm", dist.abb)
    dist.name.x <- ifelse(special, "Normal", dist.name)
    if (dist.abb == "zmlnormAlt") {
        dist.abb.x <- "lnormAlt"
        dist.name.x <- "Lognormal"
    }
    qname.x <- paste("q", dist.abb.x, sep = "")
    switch(dist.abb, lnorm = {
        x <- log(x)
    }, lnorm3 = {
        x <- log(x - param.list$threshold)
    }, zmlnorm = {
        x <- log(x[x > 0])
        censored <- censored[x > 0]
    }, zmlnormAlt = , zmnorm = {
        x <- x[x != 0]
        censored <- censored[x != 0]
    })
    param.list.x <- param.list
    if (special) {
        switch(dist.abb, lnorm = {
        }, lnorm3 = , zmlnorm = {
            param.list.x <- param.list[c("meanlog", "sdlog")]
        }, zmnorm = {
            param.list.x <- param.list[c("mean", "sd")]
        })
        names(param.list.x) <- c("mean", "sd")
    }
    if (dist.abb == "zmlnormAlt") {
        param.list.x <- param.list[c("mean", "cv")]
        names(param.list.x) <- c("mean", "cv")
    }
    dist.params.names.x <- names(param.list.x)
    ppoints.list <- ppointsCensored(x = x, censored = censored, 
        censoring.side = censoring.side, prob.method = prob.method, 
        plot.pos.con = plot.pos.con)
    q.y <- ppoints.list$Order.Statistics
    p <- ppoints.list$Cumulative.Probabilities
    cen <- ppoints.list$Censored
    q.x <- do.call(qname.x, c(list(p), param.list.x))
    if (plot.it) {
        qlab <- switch(dist.abb, lnorm = paste("Log [", x.name, 
            "]"), lnorm3 = paste("Log [", x.name, "-", format(param.list$threshold, 
            digits = digits), "]"), zmnorm = paste("Non-Zero Values of", 
            x.name), zmlnorm = paste("Log [", x.name, "> 0 ]"), 
            zmlnormAlt = paste(x.name, "> 0"), x.name)
        if (any(dist.abb == c("beta", "chisq", "f")) && param.list$ncp == 
            0) {
            x.string <- paste(dist.name, "(", paste(paste(dist.params.names[-n.dist.params], 
                signif(unlist(param.list[-n.dist.params]), digits), 
                sep = " = "), collapse = ", "), ")", sep = "")
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
            if (is.null(main)) {
                prob.method.name = switch(prob.method, `michael-schucany` = "Michael-Schucany", 
                  `hirsch-stedinger` = "Hirsch-Stedinger", `kaplan-meier` = "Kaplan-Meier", 
                  `modified kaplan-meier` = "Modified Kaplan-Meier", 
                  nelson = "Nelson")
                main <- paste(dist.name.x, " Q-Q Plot for ", 
                  qlab, ", Based on\n", prob.method.name, " Plotting Positions (Censored Data)", 
                  sep = "")
            }
        }
        else {
            if (is.null(xlab)) 
                xlab <- "Mean of Observed and Fitted Quantiles"
            if (is.null(ylab)) 
                ylab <- "Observed - Fitted Quantiles"
            if (is.null(main)) 
                main <- paste("Tukey Mean-Difference Q-Q Plot for ", 
                  qlab, " (Censored Data)\nFitted to ", x.string, 
                  " Distribution", sep = "")
        }
    }
    if (plot.it) {
        n <- length(x)
        pp.x <- ppoints(n + 2)
        emp.pct.low <- pp.x[1]
        emp.pct.up <- pp.x[n + 2]
    }
    is.finite.q.x <- is.finite(q.x)
    q.x.finite <- q.x[is.finite.q.x]
    q.y.finite <- q.y[is.finite.q.x]
    q.x.no.cen.finite <- q.x[is.finite.q.x & !cen]
    q.y.no.cen.finite <- q.y[is.finite.q.x & !cen]
    if (plot.type == "Q-Q") {
        if (plot.it) {
            if (is.null(xlim)) 
                xlim <- do.call(qname.x, c(list(c(emp.pct.low, 
                  emp.pct.up)), param.list.x))
            xlim[1] <- min(xlim[1], q.x.finite)
            xlim[2] <- max(xlim[2], q.x.finite)
            if (is.null(ylim)) {
                if (equal.axes || estimate.params) {
                  xlim[1] <- min(xlim[1], q.y.finite)
                  xlim[2] <- max(xlim[2], q.y.finite)
                  ylim <- xlim
                }
                else {
                  if (!is.null(est.param.list)) {
                    ylim <- do.call(qname, c(list(c(emp.pct.low, 
                      emp.pct.up)), est.param.list))
                    if (log.dist) 
                      ylim <- log(ylim)
                    ylim[1] <- min(ylim[1], q.y)
                    ylim[2] <- max(ylim[2], q.y)
                  }
                  else ylim <- range(q.y.finite)
                }
            }
            plot(q.x.no.cen.finite, q.y.no.cen.finite, type = "n", 
                ..., xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
                main = main)
            arg.list <- c(list(x = q.x.no.cen.finite, y = q.y.no.cen.finite, 
                method = duplicate.points.method), gen.gp.list, 
                list(col = points.col))
            do.call("points.w.dups", arg.list)
            if (add.line) 
                switch(qq.line.type, `least squares` = {
                  arg.list <- c(list(a = lm(q.y.no.cen.finite ~ 
                    q.x.no.cen.finite)), gen.gp.list, list(col = line.col, 
                    lwd = line.lwd, lty = line.lty))
                  do.call("abline", arg.list)
                }, `0-1` = {
                  arg.list <- c(list(a = 0, b = 1), gen.gp.list, 
                    list(col = line.col, lwd = line.lwd, lty = line.lty))
                  do.call("abline", arg.list)
                }, robust = {
                  arg.list <- c(list(x = q.x.no.cen.finite, y = q.y.no.cen.finite), 
                    gen.gp.list, list(col = line.col, lwd = line.lwd, 
                      lty = line.lty))
                  do.call("qqLine", arg.list)
                })
        }
        if (include.cen) {
            index <- cen & is.finite.q.x
            arg.list <- c(list(x = q.x[index], y = q.y[index], 
                method = duplicate.points.method), gen.gp.list, 
                list(pch = cen.pch, cex = cen.cex, dup.cex = cen.cex, 
                  col = cen.col))
            do.call("points.w.dups", arg.list)
            if (!all(is.finite.q.x[cen])) 
                warning(paste("Censored values corresponding to", 
                  "quantile values of Inf or -Inf not plotted.\n"))
        }
        ret.list <- c(list(x = q.x, y = q.y), ppoints.list)
    }
    else {
        q.mean <- (q.x + q.y)/2
        q.diff <- q.y - q.x
        q.mean.no.cen <- q.mean[!cen]
        q.diff.no.cen <- q.diff[!cen]
        q.mean.finite <- q.mean[is.finite.q.x]
        q.diff.finite <- q.diff[is.finite.q.x]
        if (plot.it) {
            if (is.null(ylim)) {
                rqmo2 <- diff(range(q.mean.finite))/2
                mqd <- median(q.diff.finite)
                ylim.min <- min(q.diff.finite, mqd - rqmo2)
                ylim.max <- max(q.diff.finite, mqd + rqmo2)
                ylim <- c(ylim.min, ylim.max)
            }
            if (is.null(xlim)) {
                xlim.tmp <- do.call(qname.x, c(list(c(emp.pct.low, 
                  emp.pct.up)), param.list.x))
                xlim.tmp[1] <- min(xlim.tmp[1], q.mean.finite)
                xlim.tmp[2] <- max(xlim.tmp[2], q.mean.finite)
                if (estimate.params) {
                  ylim.tmp <- xlim.tmp
                }
                else {
                  ylim.tmp <- do.call(qname, c(list(c(emp.pct.low, 
                    emp.pct.up)), param.list))
                  if (log.dist) 
                    ylim.tmp <- log(ylim.tmp)
                }
                xlim <- (xlim.tmp + ylim.tmp)/2
            }
            plot(q.mean.no.cen, q.diff.no.cen, type = "n", ..., 
                xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
                main = main)
            arg.list <- c(list(x = q.mean.no.cen, y = q.diff.no.cen, 
                method = duplicate.points.method), gen.gp.list, 
                list(col = points.col))
            do.call("points.w.dups", arg.list)
            if (add.line) {
                arg.list <- c(list(h = 0), gen.gp.list, list(col = line.col, 
                  lwd = line.lwd, lty = line.lty))
                do.call("abline", arg.list)
            }
        }
        if (include.cen) {
            index <- cen & is.finite.q.x
            arg.list <- c(list(x = q.mean[index], y = q.diff[index], 
                method = duplicate.points.method), gen.gp.list, 
                list(pch = cen.pch, cex = cen.cex, dup.cex = cen.cex, 
                  col = cen.col))
            do.call("points.w.dups", arg.list)
            if (!all(is.finite.q.x[cen])) 
                warning(paste("Censored values corresponding to", 
                  "quantile values of Inf or -Inf not plotted.\n"))
        }
        ret.list <- c(list(x = q.mean, y = q.diff), ppoints.list)
    }
    invisible(ret.list)
}
