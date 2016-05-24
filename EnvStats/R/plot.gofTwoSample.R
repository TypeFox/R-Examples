plot.gofTwoSample <-
function (x, plot.type = "Summary", captions = list(PDFs = NULL, 
    CDFs = NULL, QQ = NULL, MDQQ = NULL, Results = NULL), x.labels = list(PDFs = NULL, 
    CDFs = NULL, QQ = NULL, MDQQ = NULL), y.labels = list(PDFs = NULL, 
    CDFs = NULL, QQ = NULL, MDQQ = NULL), same.window = FALSE, 
    ask = same.window & plot.type == "All", x.points.col = "blue", 
    y.points.col = "black", points.pch = 1, jitter.points = TRUE, 
    discrete = FALSE, plot.pos.con = 0.375, x.ecdf.col = "blue", 
    y.ecdf.col = "black", x.ecdf.lwd = 3 * par("cex"), y.ecdf.lwd = 3 * 
        par("cex"), x.ecdf.lty = 1, y.ecdf.lty = 4, add.line = TRUE, 
    digits = ifelse(plot.type == "Summary", 2, .Options$digits), 
    test.result.font = 1, test.result.cex = ifelse(plot.type == 
        "Summary", 0.9, 1) * par("cex"), test.result.mar = c(0, 
        0, 3, 0) + 0.1, cex.main = ifelse(plot.type == "Summary", 
        1.2, 1.5) * par("cex"), cex.axis = ifelse(plot.type == 
        "Summary", 0.9, 1) * par("cex"), cex.lab = ifelse(plot.type == 
        "Summary", 0.9, 1) * par("cex"), main = NULL, xlab = NULL, 
    ylab = NULL, xlim = NULL, ylim = NULL, add.om.title = TRUE, 
    oma = if (plot.type == "Summary" & add.om.title) c(0, 0, 
        4, 0) else c(0, 0, 0, 0), om.title = NULL, om.font = 2, 
    om.cex.main = 1.5 * par("cex"), om.line = 0, ...) 
{
    gof.obj <- x
    plot.type <- match.arg(plot.type, c("Summary", "All", "PDFs: Observed", 
        "CDFs: Observed", "Q-Q Plot", "Tukey M-D Q-Q Plot", "Test Results"))
    check.gp.list <- checkGraphicsPars(...)
    gp.arg.list <- check.gp.list$gp.arg.list
    gen.gp.list <- check.gp.list$gen.gp.list
    data <- gof.obj$data
    x <- data[[1]]
    y <- data[[2]]
    names.data <- names(data)
    x.name <- names.data[1]
    y.name <- names.data[2]
    data.name <- gof.obj$data.name
    parent.of.data <- gof.obj$parent.of.data
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
    if (is.element(plot.type, c("Summary", "All", "PDFs: Observed"))) {
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
            else xlab <- ""
        }
        else xlab <- user.xlab
        if (is.null(user.ylab)) {
            if (!is.null(y.labels[[1]])) 
                ylab <- y.labels[[1]]
            else ylab <- "Observed Data"
        }
        else ylab <- user.ylab
        if (is.null(user.xlim)) 
            xlim <- c(0.5, 2.5)
        if (is.null(user.ylim)) 
            ylim <- range(x, y)
        dum.x.for.x <- rep(1, length(x))
        dum.x.for.y <- rep(2, length(y))
        o.mar <- par("mar")
        new.mar <- o.mar
        new.mar[2] <- max(6, o.mar[2])
        o.par <- par(mar = new.mar)
        plot(dum.x.for.x, x, type = "n", xlim = xlim, ylim = ylim, 
            axes = FALSE, xlab = xlab, ylab = ylab)
        if (jitter.points) {
            points(jitter(dum.x.for.x), x, col = x.points.col, 
                pch = points.pch)
            points(jitter(dum.x.for.y), y, col = y.points.col, 
                pch = points.pch)
        }
        else {
            points(dum.x.for.x, x, col = x.points.col, pch = points.pch)
            points(dum.x.for.y, y, col = y.points.col, pch = points.pch)
        }
        axis(1, at = 1:2, labels = c(x.name, y.name), tick = FALSE)
        axis(2)
        box()
        if (is.null(user.main)) {
            if (!is.null(captions[[1]])) 
                main <- captions[[1]]
            else main <- paste("Observed Data (Strip Plot)", 
                "\nfor", x.name, "and", y.name)
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
        par(o.par)
    }
    if (is.element(plot.type, c("Summary", "All", "CDFs: Observed"))) {
        if (plot.type == "All" & !same.window) 
            dev.new()
        if (is.null(user.xlab)) {
            if (!is.null(x.labels[[2]])) 
                xlab <- x.labels[[2]]
            else xlab <- "Order Statistics"
        }
        else xlab <- user.xlab
        if (!is.null(user.ylab)) 
            ylab <- user.ylab
        else if (!is.null(y.labels[[2]])) 
            ylab <- y.labels[[2]]
        cdfCompare(x, y, discrete = discrete, plot.pos.con = plot.pos.con, 
            x.col = x.ecdf.col, y.or.fitted.col = y.ecdf.col, 
            x.lwd = x.ecdf.lwd, y.or.fitted.lwd = y.ecdf.lwd, 
            x.lty = x.ecdf.lty, y.or.fitted.lty = y.ecdf.lty, 
            digits = digits, ..., main = "", xlab = xlab, ylab = user.ylab, 
            xlim = user.xlim, ylim = user.ylim)
        if (is.null(user.main)) {
            if (!is.null(captions[[2]])) 
                main <- captions[[2]]
            else main <- paste("Empirical CDF for ", x.name, 
                " (solid line)\nwith Empirical CDF for ", y.name, 
                " (dashed line)", sep = "", collapse = "")
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
    }
    if (is.element(plot.type, c("Summary", "All", "Q-Q Plot"))) {
        if (plot.type == "All" & !same.window) 
            dev.new()
        if (is.null(user.xlab)) {
            if (!is.null(x.labels[[3]])) 
                xlab <- x.labels[[3]]
            else xlab <- paste("Quantiles of", x.name)
        }
        else xlab <- user.xlab
        if (is.null(user.ylab)) {
            if (!is.null(y.labels[[3]])) 
                ylab <- y.labels[[3]]
            else ylab <- paste("Quantiles of", y.name)
        }
        else ylab <- user.ylab
        qqPlot(x, y, plot.pos.con = plot.pos.con, digits = digits, 
            add.line = add.line, qq.line.type = "0-1", ..., xlab = xlab, 
            ylab = ylab, main = "", xlim = user.xlim, ylim = user.ylim)
        if (is.null(user.main)) {
            if (!is.null(captions[[3]])) 
                main <- captions[[3]]
            else {
                main <- paste("Q-Q Plot of\n", y.name, "vs.", 
                  x.name)
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
        if (is.null(user.ylab)) {
            if (!is.null(y.labels[[4]])) 
                ylab <- y.labels[[4]]
            else ylab <- paste(y.name, "Quantiles -", x.name, 
                "Quantiles")
        }
        else ylab <- user.ylab
        qqPlot(x, y, plot.type = "Tukey Mean-Difference Q-Q", 
            add.line = add.line, digits = digits, ..., main = "", 
            xlab = xlab, ylab = ylab, xlim = user.xlim, ylim = user.ylim)
        if (is.null(user.main)) {
            if (!is.null(captions[[4]])) 
                main <- captions[[4]]
            else main <- paste("Tukey Mean-Difference Q-Q Plot for\n", 
                x.name, "and", y.name)
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
            else main <- paste("Results of", gof.obj$method)
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
        text(0, 1 - mf * char.ht, "Data:", adj = c(0, 1))
        text(1, 1 - mf * char.ht, paste(format(names(data.name), 
            justify = "left"), " = ", format(data.name, justify = "right"), 
            "\n", sep = "", collapse = ""), adj = c(1, 1))
        mf <- mf + 3
        if (!is.null(parent.of.data)) {
            text(0, 1 - mf * char.ht, "Data Source:", adj = 0)
            text(1, 1 - mf * char.ht, parent.of.data, adj = 1)
            mf <- mf + 2
        }
        text(0, 1 - mf * char.ht, "Sample Sizes:", adj = c(0, 
            1))
        text(1, 1 - mf * char.ht, paste(format(names(gof.obj$sample.size), 
            justify = "left"), " = ", format(format(gof.obj$sample.size, 
            digits = digits, nsmall = 0)), "\n", sep = "", collapse = ""), 
            adj = c(1, 1))
        mf <- mf + 2 + length(gof.obj$sample.size)
        text(0, 1 - mf * char.ht, "Test Statistic:", adj = 0)
        text(1, 1 - mf * char.ht, paste(names(gof.obj$statistic), 
            format(gof.obj$statistic, digits = digits, nsmall = 0), 
            sep = " = "), adj = 1)
        mf <- mf + 2
        text(0, 1 - mf * char.ht, "Test Statistic Parameters:", 
            adj = c(0, 1))
        text(1, 1 - mf * char.ht, paste(format(names(gof.obj$parameters), 
            justify = "left"), " = ", format(format(gof.obj$parameters, 
            digits = digits, nsmall = 0)), "\n", sep = "", collapse = ""), 
            adj = c(1, 1))
        mf <- mf + 2 + length(gof.obj$parameters)
        text(0, 1 - mf * char.ht, "P-value:", adj = 0)
        text(1, 1 - mf * char.ht, format(gof.obj$p.value, digits = digits, 
            nsmall = 0), adj = 1)
        par(c(o.mar, o.par))
    }
    if (plot.type == "Summary" & add.om.title) {
        if (is.null(om.title)) {
            data.name <- paste(x.name, "and", y.name)
            if (!is.null(parent.of.data)) 
                data.name <- paste(data.name, "in", parent.of.data)
            om.title <- paste("Results of ", gof.obj$method, 
                " Test for\n", data.name, sep = "")
            mtext(om.title, side = 3, line = om.line, outer = TRUE, 
                cex = om.cex.main, font = om.font)
        }
        else {
            mtext(om.title, side = 3, line = om.line, outer = TRUE, 
                cex = om.cex.main, font = om.font)
        }
    }
    invisible(gof.obj)
}
