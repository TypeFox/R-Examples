plot.gof <-
function (x, plot.type = "Summary", captions = list(PDFs = NULL, 
    CDFs = NULL, QQ = NULL, MDQQ = NULL, Results = NULL), x.labels = list(PDFs = NULL, 
    CDFs = NULL, QQ = NULL, MDQQ = NULL), y.labels = list(PDFs = NULL, 
    CDFs = NULL, QQ = NULL, MDQQ = NULL), same.window = FALSE, 
    ask = same.window & plot.type == "All", hist.col = "cyan", 
    fitted.pdf.col = "black", fitted.pdf.lwd = 3 * par("cex"), 
    fitted.pdf.lty = 1, plot.pos.con = switch(dist.abb, norm = , 
        lnorm = , lnormAlt = , lnorm3 = 0.375, evd = 0.44, 0.4), 
    ecdf.col = "cyan", fitted.cdf.col = "black", ecdf.lwd = 3 * 
        par("cex"), fitted.cdf.lwd = 3 * par("cex"), ecdf.lty = 1, 
    fitted.cdf.lty = 2, add.line = TRUE, digits = ifelse(plot.type == 
        "Summary", 2, .Options$digits), test.result.font = 1, 
    test.result.cex = ifelse(plot.type == "Summary", 0.9, 1) * 
        par("cex"), test.result.mar = c(0, 0, 3, 0) + 0.1, cex.main = ifelse(plot.type == 
        "Summary", 1.2, 1.5) * par("cex"), cex.axis = ifelse(plot.type == 
        "Summary", 0.9, 1) * par("cex"), cex.lab = ifelse(plot.type == 
        "Summary", 0.9, 1) * par("cex"), main = NULL, xlab = NULL, 
    ylab = NULL, xlim = NULL, ylim = NULL, add.om.title = TRUE, 
    oma = if (plot.type == "Summary" & add.om.title) c(0, 0, 
        2.5, 0) else c(0, 0, 0, 0), om.title = NULL, om.font = 2, 
    om.cex.main = 1.75 * par("cex"), om.line = 0.5, ...) 
{
    gof.obj <- x
    plot.type <- match.arg(plot.type, c("Summary", "All", "PDFs: Observed and Fitted", 
        "CDFs: Observed and Fitted", "Q-Q Plot", "Tukey M-D Q-Q Plot", 
        "Test Results"))
    gof.obj.orig <- gof.obj
    switch(gof.obj$dist.abb, zmnorm = gof.obj$distribution.parameters <- gof.obj$distribution.parameters[c("mean", 
        "sd", "p.zero")], zmlnorm = gof.obj$distribution.parameters <- gof.obj$distribution.parameters[c("meanlog", 
        "sdlog", "p.zero")], zmlnormAlt = gof.obj$distribution.parameters <- gof.obj$distribution.parameters[c("mean", 
        "cv", "p.zero")])
    check.gp.list <- checkGraphicsPars(...)
    gp.arg.list <- check.gp.list$gp.arg.list
    gen.gp.list <- check.gp.list$gen.gp.list
    dist.params.list <- as.list(gof.obj$distribution.parameters)
    dist.params.names <- names(gof.obj$distribution.parameters)
    names(dist.params.list) <- dist.params.names
    n.dist.params <- length(dist.params.names)
    est.method <- gof.obj$estimation.method
    dist.abb <- gof.obj$dist.abb
    dist.name <- gof.obj$distribution
    dist.type <- .Distribution.type[dist.abb]
    discrete <- !any(dist.type == c("Continuous", "Mixed"))
    data <- gof.obj$data
    data.name <- gof.obj$data.name
    data.name.string <- data.name
    parent.of.data <- gof.obj$parent.of.data
    if (!is.null(parent.of.data)) 
        data.name.string <- paste(data.name, "in", parent.of.data)
    user.main <- main
    user.xlab <- xlab
    user.ylab <- ylab
    user.xlim <- xlim
    user.ylim <- ylim
    if (!missing(captions)) {
        if (!is.list(captions)) 
            stop("The argument 'captions' must be a list")
        len <- length(captions)
        if (len < 1 | len > 5) 
            stop("The argument 'captions' must be a list with 1 to 5 components")
        if (!all(sapply(captions, length) == 1) || !all(sapply(captions, 
            is.character))) 
            stop("All components of the argument 'captions' must be character strings")
        names.vec <- names(captions)
        if (!all(names.vec %in% c("PDFs", "CDFs", "QQ", "MDQQ", 
            "Results"))) 
            stop(paste("All components of the argument 'captions'", 
                "must have names, and the name must be one of", 
                "\"PDFs\", \"CDFs\", \"QQ\", \"MDQQ\", or \"Results\""))
        if (length(unique(names.vec)) != length(names.vec)) 
            stop(paste("The names of the components for the argument 'captions'", 
                "must be unique"))
        old.captions <- captions
        captions <- list(PDFs = NULL, CDFs = NULL, QQ = NULL, 
            MDQQ = NULL, Results = NULL)
        captions[names.vec] <- old.captions
    }
    if (!missing(x.labels)) {
        if (!is.list(x.labels)) 
            stop("The argument 'x.labels' must be a list")
        len <- length(x.labels)
        if (len < 1 | len > 4) 
            stop("The argument 'x.labels' must be a list with 1 to 4 components")
        if (!all(sapply(x.labels, length) == 1) || !all(sapply(x.labels, 
            is.character))) 
            stop("All components of the argument 'x.labels' must be character strings")
        names.vec <- names(x.labels)
        if (!all(names.vec %in% c("PDFs", "CDFs", "QQ", "MDQQ"))) 
            stop(paste("All components of the argument 'x.labels'", 
                "must have names, and the name must be one of", 
                "\"PDFs\", \"CDFs\", \"QQ\", or \"MDQQ\""))
        if (length(unique(names.vec)) != length(names.vec)) 
            stop(paste("The names of the components for the argument 'x.labels'", 
                "must be unique"))
        old.x.labels <- x.labels
        x.labels <- list(PDFs = NULL, CDFs = NULL, QQ = NULL, 
            MDQQ = NULL)
        x.labels[names.vec] <- old.x.labels
    }
    if (!missing(y.labels)) {
        if (!is.list(y.labels)) 
            stop("The argument 'y.labels' must be a list")
        len <- length(y.labels)
        if (len < 1 | len > 4) 
            stop("The argument 'y.labels' must be a list with 1 to 4 components")
        if (!all(sapply(y.labels, length) == 1) || !all(sapply(y.labels, 
            is.character))) 
            stop("All components of the argument 'y.labels' must be character strings")
        names.vec <- names(y.labels)
        if (!all(names.vec %in% c("PDFs", "CDFs", "QQ", "MDQQ"))) 
            stop(paste("All components of the argument 'y.labels'", 
                "must have names, and the name must be one of", 
                "\"PDFs\", \"CDFs\", \"QQ\", or \"MDQQ\""))
        if (length(unique(names.vec)) != length(names.vec)) 
            stop(paste("The names of the components for the argument 'y.labels'", 
                "must be unique"))
        old.y.labels <- y.labels
        y.labels <- list(PDFs = NULL, CDFs = NULL, QQ = NULL, 
            MDQQ = NULL)
        y.labels[names.vec] <- old.y.labels
    }
    if (is.element(plot.type, c("Summary", "All", "PDFs: Observed and Fitted"))) {
        if (plot.type == "All" & same.window) {
            devAskNewPage(ask = ask)
        }
        else if (plot.type == "Summary") {
            o.par1 <- par(c("cex", "mex", "mgp"))
            o.par2 <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1) + 
                0.1, oma = oma)
            par(cex = 0.8 * o.par1$cex, mex = 0.8 * o.par1$mex, 
                mgp = c(2.75, 0.5, 0))
            on.exit(par(c(o.par1, o.par2)))
        }
        if (is.null(user.xlab)) {
            if (!is.null(x.labels[[1]])) 
                xlab <- x.labels[[1]]
            else xlab <- data.name
        }
        else xlab <- user.xlab
        if (is.null(user.ylab)) {
            if (!is.null(y.labels[[1]])) 
                ylab <- y.labels[[1]]
            else ylab <- "Relative Frequency"
        }
        else ylab <- user.ylab
        if (!discrete) {
            if (!is.null(gof.obj$cut.points) && all(is.finite(range(gof.obj$cut.points)))) {
                hist.list <- hist(data, breaks = gof.obj$cut.points, 
                  plot = FALSE)
                if (is.null(user.xlim)) 
                  xlim <- range(hist.list$breaks)
                pdf.list <- do.call("pdfPlot", list(distribution = dist.abb, 
                  param.list = dist.params.list, plot.it = FALSE, 
                  ..., xlim = xlim))
                if (is.null(user.ylim)) 
                  ylim <- c(0, max(hist.list$density, pdf.list$Probability.Densities))
                hist(data, probability = TRUE, breaks = gof.obj$cut.points, 
                  col = hist.col, main = "", ..., cex.axis = cex.axis, 
                  cex.lab = cex.lab, xlim = xlim, ylim = ylim, 
                  xlab = xlab, ylab = ylab)
            }
            else {
                data.no.0 <- data[data != 0]
                if (dist.type == "Mixed" && all(data.no.0 > 0)) {
                  upper.breaks <- pretty(data.no.0, n = log(length(data.no.0), 
                    base = 2) + 1)
                  if (min(data.no.0) <= upper.breaks[1]) 
                    upper.breaks[1] <- upper.breaks[1] - 1e+08 * 
                      .Machine$double.eps
                  if (upper.breaks[1] > 0) {
                    upper.0.break <- min(0.5, upper.breaks[1] - 
                      1e+08 * .Machine$double.eps)
                    breaks <- c(-upper.0.break, upper.0.break, 
                      upper.breaks)
                  }
                  else breaks <- upper.breaks
                }
                else {
                  breaks <- pretty(data, n = log(length(data), 
                    base = 2) + 1)
                }
                hist.list <- hist(data, breaks = breaks, plot = FALSE)
                if (is.null(user.xlim)) 
                  xlim <- range(hist.list$breaks)
                pdf.list <- do.call("pdfPlot", list(distribution = dist.abb, 
                  param.list = dist.params.list, plot.it = FALSE, 
                  ..., xlim = xlim))
                if (is.null(user.ylim)) 
                  ylim <- c(0, max(hist.list$density, pdf.list$Probability.Densities))
                hist(data, breaks = breaks, probability = TRUE, 
                  col = hist.col, main = "", ..., cex.axis = cex.axis, 
                  cex.lab = cex.lab, xlim = xlim, ylim = ylim, 
                  xlab = xlab, ylab = ylab)
            }
        }
        else {
            y <- tabulate(data - min(data) + 1)/length(data)
            pdf.list <- do.call("pdfPlot", list(distribution = dist.abb, 
                param.list = dist.params.list, plot.it = FALSE, 
                ..., xlim = xlim))
            if (is.null(user.ylim)) 
                ylim <- c(0, max(y, pdf.list$Probability.Densities))
            x <- min(data):max(data)
            nx <- length(x)
            con <- 0.4 + (0.1 * (nx - 2))/nx
            xleft <- x - con
            xright <- x + con
            ybottom <- rep(0, nx)
            if (is.null(user.xlim)) 
                xlim <- c(min(xleft), max(xright))
            plot(x, y, type = "n", xaxt = "n", bty = "n", ..., 
                cex.axis = cex.axis, cex.lab = cex.lab, xlim = xlim, 
                ylim = ylim, xlab = xlab, ylab = ylab)
            rect(xleft = xleft, ybottom = ybottom, xright = xright, 
                ytop = y, col = hist.col, border = fitted.pdf.col, 
                ...)
            axis(1)
        }
        arg.list <- c(list(distribution = dist.abb, param.list = dist.params.list, 
            add = TRUE, pdf.col = fitted.pdf.col, pdf.lwd = fitted.pdf.lwd, 
            pdf.lty = fitted.pdf.lty), gen.gp.list)
        do.call("pdfPlot", arg.list)
        if (is.null(user.main)) {
            if (!is.null(captions[[1]])) 
                main <- captions[[1]]
            else main <- paste("Histogram for ", data.name, " with\nFitted ", 
                gof.obj$distribution, " Distribution", sep = "")
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
    }
    if (is.element(plot.type, c("Summary", "All", "CDFs: Observed and Fitted"))) {
        if (plot.type == "All" & !same.window) 
            dev.new()
        if (is.null(user.xlab)) {
            if (!is.null(x.labels[[2]])) 
                xlab <- x.labels[[2]]
            else {
                if (is.null(est.method)) {
                  string <- dist.name
                }
                else {
                  if (any(dist.abb == c("beta", "chisq", "f", 
                    "t"))) {
                    if (dist.params.list$ncp > 0) 
                      string <- paste("Non-central ", dist.name, 
                        "(", paste(paste(dist.params.names, signif(unlist(dist.params.list), 
                          digits), sep = "="), collapse = ", "), 
                        ")", sep = "")
                    else {
                      string <- paste(dist.name, "(", paste(paste(dist.params.names[-n.dist.params], 
                        signif(unlist(dist.params.list[-n.dist.params]), 
                          digits), sep = "="), collapse = ", "), 
                        ")", sep = "")
                    }
                  }
                  else {
                    string <- paste(dist.name, "(", paste(paste(dist.params.names, 
                      signif(unlist(dist.params.list), digits), 
                      sep = "="), collapse = ", "), ")", sep = "")
                  }
                }
                xlab <- paste("Order Statistics for ", data.name, 
                  " and\n", string, " Distribution", sep = "")
            }
        }
        else xlab <- user.xlab
        if (!is.null(user.ylab)) 
            ylab <- user.ylab
        else if (!is.null(y.labels[[2]])) 
            ylab <- y.labels[[2]]
        cdfCompare(data, plot.pos.con = plot.pos.con, distribution = dist.abb, 
            param.list = dist.params.list, estimate.params = FALSE, 
            x.col = ecdf.col, y.or.fitted.col = fitted.cdf.col, 
            x.lwd = ecdf.lwd, y.or.fitted.lwd = fitted.cdf.lwd, 
            x.lty = ecdf.lty, y.or.fitted.lty = fitted.cdf.lty, 
            digits = digits, ..., cex.axis = cex.axis, cex.lab = cex.lab, 
            main = "", xlab = xlab, ylab = ylab, xlim = user.xlim, 
            ylim = user.ylim)
        if (is.null(user.main)) {
            if (!is.null(captions[[2]])) 
                main <- captions[[2]]
            else main <- paste("Empirical CDF for ", data.name, 
                " (solid line)\nwith Fitted ", gof.obj$distribution, 
                " CDF (dashed line)", sep = "", collapse = "")
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
    }
    if (is.element(plot.type, c("Summary", "All", "Q-Q Plot"))) {
        if (plot.type == "All" & !same.window) 
            dev.new()
        if (is.null(user.xlab) & !is.null(x.labels[[3]])) 
            xlab <- x.labels[[3]]
        else xlab <- user.xlab
        if (is.null(user.ylab)) {
            if (!is.null(y.labels[[3]])) 
                ylab <- y.labels[[3]]
            else {
                qlab <- switch(dist.abb, lnorm = paste("Log[", 
                  data.name, "]", sep = ""), lnorm3 = paste("Log[", 
                  data.name, "-", format(dist.params.list[["threshold"]], 
                    digits = digits), "]", sep = ""), data.name)
                ylab <- paste("Quantiles of", qlab)
            }
        }
        else ylab <- user.ylab
        qqPlot(data, distribution = dist.abb, param.list = dist.params.list, 
            plot.pos.con = plot.pos.con, digits = digits, add.line = add.line, 
            qq.line.type = "0-1", ..., cex.axis = cex.axis, cex.lab = cex.lab, 
            xlab = xlab, ylab = ylab, main = "", xlim = user.xlim, 
            ylim = user.ylim)
        if (is.null(user.main)) {
            if (!is.null(captions[[3]])) 
                main <- captions[[3]]
            else {
                main <- paste("Q-Q Plot for ", data.name, " Fitted to\n", 
                  gof.obj$distribution, " Distribution", sep = "")
                if (add.line) 
                  main <- paste(main, ", with 0-1 Line", sep = "")
            }
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
    }
    if (is.element(plot.type, c("All", "Tukey M-D Q-Q Plot"))) {
        if (plot.type == "All" & !same.window) 
            dev.new()
        if (is.null(user.xlab) & !is.null(x.labels[[4]])) 
            xlab <- x.labels[[4]]
        else xlab <- user.xlab
        if (is.null(user.ylab) & !is.null(y.labels[[4]])) 
            ylab <- y.labels[[4]]
        else ylab <- user.ylab
        qqPlot(data, distribution = dist.abb, param.list = dist.params.list, 
            plot.type = "Tukey Mean-Difference Q-Q", add.line = add.line, 
            digits = digits, ..., cex.axis = cex.axis, cex.lab = cex.lab, 
            main = "", xlab = xlab, ylab = ylab, xlim = user.xlim, 
            ylim = user.ylim)
        if (is.null(user.main)) {
            if (!is.null(captions[[4]])) 
                main <- captions[[4]]
            else main <- paste("Tukey Mean-Difference Q-Q Plot\nfor ", 
                data.name, " Fitted to ", gof.obj$distribution, 
                " Distribution", sep = "")
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
    }
    if (is.element(plot.type, c("Summary", "All", "Test Results"))) {
        if (plot.type == "All" & !same.window) 
            dev.new()
        par(usr = c(0, 1, 0, 1))
        o.mar <- par(mar = test.result.mar)
        plot(0:1, 0:1, type = "n", axes = FALSE, main = "")
        if (is.null(user.main)) {
            if (!is.null(captions[[5]])) 
                main <- captions[[5]]
            else main <- gof.obj$method
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
        o.par <- par(cex = test.result.cex, font = test.result.font)
        char.ht <- par("cxy")[2] * test.result.cex
        text(0, 1, "Hypothesized\nDistribution:", adj = c(0, 
            1))
        text(1, 1 - char.ht, gof.obj$distribution, adj = c(1, 
            1))
        mf <- 3
        if (gof.obj$n.param.est > 0) {
            text(0, 1 - mf * char.ht, "Estimated Parameters:", 
                adj = c(0, 1))
            text(1, 1 - mf * char.ht, paste(format(names(gof.obj$distribution.parameters), 
                justify = "left"), " = ", format(format(gof.obj$distribution.parameters, 
                digits = digits, nsmall = 0)), "\n", sep = "", 
                collapse = ""), adj = c(1, 1))
            mf <- mf + 2 + length(gof.obj$distribution.parameters) - 
                1
        }
        text(0, 1 - mf * char.ht, "Data:", adj = 0)
        text(1, 1 - mf * char.ht, gof.obj$data.name, adj = 1)
        mf <- mf + 2
        if (!is.null(gof.obj$parent.of.data)) {
            text(0, 1 - mf * char.ht, "Data Source:", adj = 0)
            text(1, 1 - mf * char.ht, gof.obj$parent.of.data, 
                adj = 1)
            mf <- mf + 2
        }
        text(0, 1 - mf * char.ht, "Sample Size:", adj = 0)
        text(1, 1 - mf * char.ht, gof.obj$sample.size, adj = 1)
        mf <- mf + 2
        text(0, 1 - mf * char.ht, "Test Statistic:", adj = 0)
        text(1, 1 - mf * char.ht, paste(names(gof.obj$statistic), 
            format(gof.obj$statistic, digits = digits, nsmall = 0), 
            sep = " = "), adj = 1)
        mf <- mf + 2
        n.params <- length(gof.obj$parameters)
        if (n.params > 0) {
            text(0, 1 - mf * char.ht, "Test Statistic Parameter:", 
                adj = 0)
            text(1, 1 - mf * char.ht, paste(format(names(gof.obj$parameters), 
                justify = "left"), format(gof.obj$parameters, 
                digits = digits, nsmall = 0), sep = " = "), adj = 1)
            mf <- mf + 2
        }
        text(0, 1 - mf * char.ht, "P-value:", adj = 0)
        text(1, 1 - mf * char.ht, format(gof.obj$p.value, digits = digits, 
            nsmall = 0), adj = 1)
        par(c(o.mar, o.par))
    }
    if (plot.type == "Summary" & add.om.title) {
        if (is.null(om.title)) 
            om.title <- paste("Goodness-of-Fit Results for", 
                data.name.string)
        mtext(om.title, side = 3, line = om.line, outer = TRUE, 
            cex = om.cex.main, font = om.font)
    }
    invisible(gof.obj.orig)
}
