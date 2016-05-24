plot.gofCensored <-
function (x, plot.type = "Summary", captions = list(PDFs = NULL, 
    CDFs = NULL, QQ = NULL, MDQQ = NULL, Results = NULL), x.labels = list(PDFs = NULL, 
    CDFs = NULL, QQ = NULL, MDQQ = NULL), y.labels = list(PDFs = NULL, 
    CDFs = NULL, QQ = NULL, MDQQ = NULL), same.window = FALSE, 
    ask = same.window & plot.type == "All", hist.col = "cyan", 
    fitted.pdf.col = "black", fitted.pdf.lwd = 3 * par("cex"), 
    fitted.pdf.lty = 1, prob.method = "michael-schucany", plot.pos.con = 0.375, 
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
        4, 0) else c(0, 0, 0, 0), om.title = NULL, om.font = 2, 
    om.cex.main = 1.5 * par("cex"), om.line = 0, ...) 
{
    gofCensored.obj <- x
    plot.type <- match.arg(plot.type, c("Summary", "All", "PDFs: Observed and Fitted", 
        "CDFs: Observed and Fitted", "Q-Q Plot", "Tukey M-D Q-Q Plot", 
        "Test Results"))
    check.gp.list <- checkGraphicsPars(...)
    gp.arg.list <- check.gp.list$gp.arg.list
    gen.gp.list <- check.gp.list$gen.gp.list
    dist.params.list <- as.list(gofCensored.obj$distribution.parameters)
    dist.params.names <- names(gofCensored.obj$distribution.parameters)
    names(dist.params.list) <- dist.params.names
    n.dist.params <- length(dist.params.names)
    dist.abb <- gofCensored.obj$dist.abb
    dist.name <- gofCensored.obj$distribution
    dist.type <- .Distribution.type[dist.abb]
    discrete <- !any(dist.type == c("Continuous", "Mixed"))
    data <- gofCensored.obj$data
    data.name <- gofCensored.obj$data.name
    data.name.string <- data.name
    parent.of.data <- gofCensored.obj$parent.of.data
    if (!is.null(parent.of.data)) 
        data.name.string <- paste(data.name, "in", parent.of.data)
    censored <- gofCensored.obj$censored
    censoring.side <- gofCensored.obj$censoring.side
    censoring.levels <- gofCensored.obj$censoring.levels
    data.no.cen <- data[!censored]
    data.cen <- data[censored]
    prob.method <- match.arg(prob.method, c("michael-schucany", 
        "hirsch-stedinger", "kaplan-meier", "modified kaplan-meier", 
        "nelson"))
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
        if (length(censoring.levels) > 1 && ((censoring.side == 
            "left" && any(data.no.cen < max(censoring.levels)))) || 
            ((censoring.side == "right" && any(data.no.cen > 
                min(censoring.levels))))) {
            warning(paste("Cannot construct histogram for", "multiply censored data when", 
                "complete observations are between", "censoring levels"))
        }
        else {
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
            if (is.null(user.xlim)) {
                qname <- paste("q", dist.abb, sep = "")
                xlim <- do.call(qname, c(list(p = c(0.001, 0.999)), 
                  dist.params.list))
                if (censoring.side == "left") {
                  xlim[1] <- min(xlim[1], min(data.cen) - 1e+08 * 
                    .Machine$double.eps)
                  xlim[2] <- max(xlim[2], max(data.no.cen))
                }
                else {
                  xlim[1] <- min(xlim[1], min(data.no.cen))
                  xlim[2] <- max(xlim[2], max(data.cen) + 1e+08 * 
                    .Machine$double.eps)
                }
            }
            pdf.list <- do.call("pdfPlot", list(distribution = dist.abb, 
                param.list = dist.params.list, plot.it = FALSE, 
                xlim = xlim, ...))
            if (!discrete) {
                if (censoring.side == "left") {
                  mcl <- max(censoring.levels)
                  data.for.breaks <- c(mcl, data.no.cen)
                  breaks <- pretty(range(data.for.breaks), n = nclass.Sturges(data.for.breaks))
                  breaks[1] <- mcl
                  breaks <- c(xlim[1], breaks)
                  data.cen.for.hist <- data.cen
                  data.cen.for.hist[data.cen.for.hist == mcl] <- mcl - 
                    1e-07 * stats::median(diff(breaks)) - 1e+08 * 
                    .Machine$double.eps
                }
                else {
                  mcl <- min(censoring.levels)
                  data.for.breaks <- c(data.no.cen, mcl)
                  breaks <- pretty(range(data.for.breaks), n = nclass.Sturges(data.for.breaks))
                  breaks[length(breaks)] <- mcl
                  breaks <- c(breaks, xlim[2])
                  data.cen.for.hist <- data.cen
                  data.cen.for.hist[data.cen.for.hist == mcl] <- mcl + 
                    1e-07 * stats::median(diff(breaks)) + 1e+08 * 
                    .Machine$double.eps
                }
                data.for.hist <- c(data.cen.for.hist, data.no.cen)
                hist.list <- hist(data.for.hist, breaks = breaks, 
                  plot = FALSE)
                if (is.null(user.ylim)) {
                  ylim <- range(pretty(c(0, max(hist.list$density, 
                    pdf.list$Probability.Densities))))
                  ylim[1] <- 0
                }
                hist(data.for.hist, breaks = breaks, right = censoring.side == 
                  "right", probability = TRUE, col = hist.col, 
                  main = "", cex.axis = cex.axis, cex.lab = cex.lab, 
                  xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, 
                  ...)
            }
            else {
                n <- length(data)
                props.no.cen <- tabulate(data.no.cen - min(data.no.cen) + 
                  1)/n
                props.cen <- length(data.cen)/n
                if (censoring.side == "left") {
                  x <- c(max(data.cen), min(data.no.cen):max(data.no.cen))
                  y <- c(props.cen, props.no.cen)
                }
                else {
                  x <- c(min(data.no.cen):max(data.no.cen), min(data.cen))
                  y <- c(props.no.cen, props.cen)
                }
                if (is.null(user.ylim)) 
                  ylim <- c(0, max(y, pdf.list$Probability.Densities))
                nx <- length(x)
                con <- 0.4 + (0.1 * (nx - 2))/nx
                xleft <- x - con
                xright <- x + con
                ybottom <- rep(0, nx)
                if (is.null(user.xlim)) 
                  xlim <- c(min(xleft), max(xright))
                plot(x, y, type = "n", xaxt = "n", bty = "n", 
                  cex.axis = cex.axis, cex.lab = cex.lab, xlim = xlim, 
                  ylim = ylim, xlab = xlab, ylab = ylab, ...)
                rect(xleft = xleft, ybottom = ybottom, xright = xright, 
                  ytop = y, col = hist.col, border = fitted.pdf.col, 
                  ...)
                axis(1, cex.axis = cex.axis, cex.lab = cex.lab)
            }
            arg.list <- c(list(distribution = dist.abb, param.list = dist.params.list, 
                add = TRUE, pdf.col = fitted.pdf.col, pdf.lwd = fitted.pdf.lwd, 
                pdf.lty = fitted.pdf.lty), gen.gp.list)
            do.call("pdfPlot", arg.list)
            if (is.null(user.main)) {
                if (!is.null(captions[[1]])) 
                  main <- captions[[1]]
                else main <- paste("Histogram for ", data.name, 
                  "\nwith Fitted ", gofCensored.obj$distribution, 
                  " Distribution", sep = "")
            }
            arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
                list(main = main))
            do.call("title", arg.list)
        }
    }
    if (is.element(plot.type, c("Summary", "All", "CDFs: Observed and Fitted"))) {
        if (plot.type == "All" & !same.window) 
            dev.new()
        if (is.null(user.xlab)) {
            if (!is.null(x.labels[[2]])) 
                xlab <- x.labels[[2]]
            else {
                if (any(dist.abb == c("beta", "chisq", "f", "t"))) {
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
                xlab <- paste("Order Statistics for ", data.name, 
                  " and\n", string, " Distribution", sep = "")
            }
        }
        else xlab <- user.xlab
        if (!is.null(user.ylab)) 
            ylab <- user.ylab
        else if (!is.null(y.labels[[2]])) 
            ylab <- y.labels[[2]]
        cdfCompareCensored(x = data, censored = censored, censoring.side = censoring.side, 
            prob.method = prob.method, plot.pos.con = plot.pos.con, 
            distribution = dist.abb, param.list = dist.params.list, 
            estimate.params = FALSE, x.col = ecdf.col, y.or.fitted.col = fitted.cdf.col, 
            x.lwd = ecdf.lwd, y.or.fitted.lwd = fitted.cdf.lwd, 
            x.lty = ecdf.lty, y.or.fitted.lty = fitted.cdf.lty, 
            digits = digits, cex.axis = cex.axis, cex.lab = cex.lab, 
            main = "", xlab = xlab, ylab = user.ylab, xlim = user.xlim, 
            ylim = user.ylim, ...)
        if (is.null(user.main)) {
            if (!is.null(captions[[2]])) 
                main <- captions[[2]]
            else main <- paste("Empirical CDF for ", data.name, 
                " (solid line)\nwith Fitted ", gofCensored.obj$distribution, 
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
        qqPlotCensored(x = data, censored = censored, censoring.side = censoring.side, 
            prob.method = prob.method, plot.pos.con = plot.pos.con, 
            distribution = dist.abb, param.list = dist.params.list, 
            digits = digits, add.line = add.line, qq.line.type = "0-1", 
            cex.axis = cex.axis, cex.lab = cex.lab, main = "", 
            xlab = user.xlab, ylab = ylab, xlim = user.xlim, 
            ylim = user.ylim, ...)
        if (is.null(user.main)) {
            if (!is.null(captions[[3]])) 
                main <- captions[[3]]
            else {
                main <- paste("Q-Q Plot for", data.name, "\nFitted to", 
                  gofCensored.obj$distribution, "Distribution")
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
        qqPlotCensored(x = data, censored = censored, censoring.side = censoring.side, 
            prob.method = prob.method, plot.pos.con = plot.pos.con, 
            distribution = dist.abb, param.list = dist.params.list, 
            plot.type = "Tukey Mean-Difference Q-Q", add.line = add.line, 
            digits = digits, ..., cex.axis = cex.axis, cex.lab = cex.lab, 
            main = "", xlab = user.xlab, ylab = user.ylab, xlim = user.xlim, 
            ylim = user.ylim)
        if (is.null(user.main)) {
            if (!is.null(captions[[4]])) 
                main <- captions[[4]]
            else main <- paste("Tukey Mean-Difference Q-Q Plot\nfor ", 
                data.name, " Fitted to ", gofCensored.obj$distribution, 
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
        method <- gofCensored.obj$method
        if (is.null(user.main)) {
            if (!is.null(captions[[5]])) 
                main <- captions[[5]]
            else {
                strings <- unlist(strsplit(method, "\n", fixed = TRUE))
                if (length(strings) > 1) {
                  string1 <- strings[1]
                  string2 <- strings[2]
                  string2 <- substring(string2, 34, nchar(string2))
                  main <- paste("Results of ", string1, "\n", 
                    string2, sep = "")
                }
                else {
                  main <- paste("Results of", method)
                }
            }
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
        o.par <- par(cex = test.result.cex, font = test.result.font)
        char.ht <- par("cxy")[2] * test.result.cex
        text(0, 1, "Hypothesized\nDistribution:", adj = c(0, 
            1))
        text(1, 1 - char.ht, gofCensored.obj$distribution, adj = c(1, 
            1))
        mf <- 3
        if (gofCensored.obj$n.param.est > 0) {
            text(0, 1 - mf * char.ht, "Estimated Parameters:", 
                adj = c(0, 1))
            text(1, 1 - mf * char.ht, paste(format(names(gofCensored.obj$distribution.parameters), 
                justify = "left"), " = ", format(gofCensored.obj$distribution.parameters, 
                digits = digits, nsmall = 0), "\n", sep = "", 
                collapse = ""), adj = c(1, 1))
            mf <- mf + 2 + length(gofCensored.obj$distribution.parameters) - 
                1
        }
        text(0, 1 - mf * char.ht, "Data:", adj = 0)
        text(1, 1 - mf * char.ht, data.name, adj = 1)
        mf <- mf + 2
        text(0, 1 - mf * char.ht, "Sample Size:", adj = 0)
        text(1, 1 - mf * char.ht, gofCensored.obj$sample.size, 
            adj = 1)
        mf <- mf + 2
        text(0, 1 - mf * char.ht, "Test Statistic:", adj = 0)
        text(1, 1 - mf * char.ht, paste(names(gofCensored.obj$statistic), 
            format(gofCensored.obj$statistic, digits = digits, 
                nsmall = 0), sep = " = "), adj = 1)
        mf <- mf + 2
        text(0, 1 - mf * char.ht, "Test Statistic Parmeters:", 
            adj = c(0, 1))
        text(1, 1 - mf * char.ht, paste(format(names(gofCensored.obj$parameters), 
            justify = "left"), " = ", format(gofCensored.obj$parameters, 
            digits = digits, nsmall = 0), "\n", sep = "", collapse = ""), 
            adj = c(1, 1))
        mf <- mf + 2 + length(gofCensored.obj$parameters) - 1
        text(0, 1 - mf * char.ht, "P-value:", adj = 0)
        text(1, 1 - mf * char.ht, format(gofCensored.obj$p.value, 
            digits = digits, nsmall = 0), adj = 1)
        par(c(o.mar, o.par))
    }
    if (plot.type == "Summary" & add.om.title) {
        if (is.null(om.title)) {
            data.name <- gofCensored.obj$data.name
            parent.of.data <- gofCensored.obj$parent.of.data
            if (!is.null(parent.of.data)) 
                data.name <- paste(data.name, "in", parent.of.data)
            method <- gofCensored.obj$method
            strings <- unlist(strsplit(method, "\n", fixed = TRUE))
            if (length(strings) > 1) {
                string1 <- strings[1]
                string2 <- strings[2]
                string2 <- substring(string2, 34, nchar(string2))
                om.title <- paste("Results of", string1, string2, 
                  "Test\nfor", data.name)
            }
            else {
                om.title <- paste("Results of", method, "Test\nfor", 
                  data.name)
            }
        }
        mtext(om.title, side = 3, line = om.line, outer = TRUE, 
            cex = om.cex.main, font = om.font)
    }
    invisible(gofCensored.obj)
}
