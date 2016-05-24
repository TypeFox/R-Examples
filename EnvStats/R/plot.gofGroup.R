plot.gofGroup <-
function (x, plot.type = "Summary", captions = list(QQ = NULL, 
    MDQQ = NULL, ScoresQQ = NULL, ScoresMDQQ = NULL, Results = NULL), 
    x.labels = list(QQ = NULL, MDQQ = NULL, ScoresQQ = NULL, 
        ScoresMDQQ = NULL), y.labels = list(QQ = NULL, MDQQ = NULL, 
        ScoresQQ = NULL, ScoresMDQQ = NULL), same.window = FALSE, 
    ask = same.window & plot.type == "All", add.line = TRUE, 
    digits = ifelse(plot.type == "Summary", 2, .Options$digits), 
    test.result.font = 1, test.result.cex = ifelse(plot.type == 
        "Summary", 0.9, 1) * par("cex"), test.result.mar = c(0, 
        0, 3, 0) + 0.1, individual.p.values = FALSE, cex.main = ifelse(plot.type == 
        "Summary", 1.2, 1.5) * par("cex"), cex.axis = ifelse(plot.type == 
        "Summary", 0.9, 1) * par("cex"), cex.lab = ifelse(plot.type == 
        "Summary", 0.9, 1) * par("cex"), main = NULL, xlab = NULL, 
    ylab = NULL, xlim = NULL, ylim = NULL, add.om.title = TRUE, 
    oma = if (plot.type == "Summary" & add.om.title) c(0, 0, 
        5, 0) else c(0, 0, 0, 0), om.title = NULL, om.font = 2, 
    om.cex.main = 1.5 * par("cex"), om.line = 1, ...) 
{
    gof.obj <- x
    plot.type <- match.arg(plot.type, c("Summary", "All", "Q-Q Plot", 
        "Tukey M-D Q-Q Plot", "Scores Q-Q Plot", "Scores Tukey M-D Q-Q Plot", 
        "Test Results"))
    check.gp.list <- checkGraphicsPars(...)
    gp.arg.list <- check.gp.list$gp.arg.list
    gen.gp.list <- check.gp.list$gen.gp.list
    data.name <- gof.obj$data.name
    data.name.string <- data.name
    parent.of.data <- gof.obj$parent.of.data
    if (!is.null(parent.of.data)) 
        data.name.string <- paste(data.name, "in", parent.of.data)
    grouping.variable <- gof.obj$grouping.variable
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
        if (!all(names.vec %in% c("QQ", "MDQQ", "ScoresQQ", "ScoresMDQQ", 
            "Results"))) 
            stop(paste("All components of the argument 'captions'", 
                "must have names, and the name must be one of", 
                "\"QQ\", \"MDQQ\", \"ScoresQQ\", \"ScoresMDQQ\", or \"Results\""))
        if (length(unique(names.vec)) != length(names.vec)) 
            stop(paste("The names of the components for the argument 'captions'", 
                "must be unique"))
        old.captions <- captions
        captions <- list(QQ = NULL, MDQQ = NULL, ScoresQQ = NULL, 
            ScoresMDQQ = NULL, Results = NULL)
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
        if (!all(names.vec %in% c("QQ", "MDQQ", "ScoresQQ", "ScoresMDQQ"))) 
            stop(paste("All components of the argument 'x.labels'", 
                "must have names, and the name must be one of", 
                "\"QQ\", \"MDQQ\", \"ScoresQQ\", or \"ScoresMDQQ\""))
        if (length(unique(names.vec)) != length(names.vec)) 
            stop(paste("The names of the components for the argument 'x.labels'", 
                "must be unique"))
        old.x.labels <- x.labels
        x.labels <- list(QQ = NULL, MDQQ = NULL, ScoresQQ = NULL, 
            ScoresMDQQ = NULL)
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
        if (!all(names.vec %in% c("QQ", "MDQQ", "ScoresQQ", "ScoresMDQQ"))) 
            stop(paste("All components of the argument 'y.labels'", 
                "must have names, and the name must be one of", 
                "\"QQ\", \"MDQQ\", \"ScoresQQ\", or \"ScoresMDQQ\""))
        if (length(unique(names.vec)) != length(names.vec)) 
            stop(paste("The names of the components for the argument 'y.labels'", 
                "must be unique"))
        old.y.labels <- y.labels
        y.labels <- list(QQ = NULL, MDQQ = NULL, ScoresQQ = NULL, 
            ScoresMDQQ = NULL)
        y.labels[names.vec] <- old.y.labels
    }
    if (is.element(plot.type, c("Summary", "All", "Q-Q Plot"))) {
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
        if (is.null(user.xlab) & !is.null(x.labels[[1]])) 
            xlab <- x.labels[[1]]
        else xlab <- user.xlab
        if (is.null(user.ylab)) {
            if (!is.null(y.labels[[1]])) 
                ylab <- y.labels[[1]]
            else {
                ylab <- paste("Quantiles of P-values")
            }
        }
        else ylab <- user.ylab
        qqPlot(gof.obj$p.value[1:gof.obj$n.groups], distribution = "unif", 
            param.list = list(min = 0, max = 1), digits = digits, 
            add.line = add.line, qq.line.type = "0-1", ..., cex.axis = cex.axis, 
            cex.lab = cex.lab, xlab = user.xlab, ylab = ylab, 
            main = "", xlim = user.xlim, ylim = user.ylim)
        if (is.null(user.main)) {
            if (!is.null(captions[[1]])) 
                main <- captions[[1]]
            else {
                main <- paste("Quantile-Quantile Plot of P-values for\n", 
                  paste(data.name, "Grouped by", grouping.variable))
            }
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
    }
    if (is.element(plot.type, c("All", "Tukey M-D Q-Q Plot"))) {
        if (plot.type == "All" & !same.window) 
            dev.new()
        if (is.null(user.xlab) & !is.null(x.labels[[2]])) 
            xlab <- x.labels[[2]]
        else xlab <- user.xlab
        if (is.null(user.ylab) & !is.null(y.labels[[2]])) 
            ylab <- y.labels[[2]]
        else ylab <- user.ylab
        qqPlot(gof.obj$p.value[1:gof.obj$n.groups], distribution = "unif", 
            param.list = list(min = 0, max = 1), plot.type = "Tukey Mean-Difference Q-Q", 
            digits = digits, add.line = add.line, ..., cex.axis = cex.axis, 
            cex.lab = cex.lab, xlab = user.xlab, ylab = ylab, 
            main = "", xlim = user.xlim, ylim = user.ylim)
        if (is.null(user.main)) {
            if (!is.null(captions[[2]])) 
                main <- captions[[2]]
            else main <- paste("Tukey Mean-Difference Q-Q Plot of P-values for\n", 
                paste(data.name, "Grouped by", grouping.variable))
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
    }
    if (is.element(plot.type, c("Summary", "All", "Scores Q-Q Plot"))) {
        method <- gof.obj$method
        n.char <- nchar(method)
        score.type <- ifelse(substr(method, n.char - 13, n.char - 
            1) == "Normal Scores", "Normal", "Chi-square")
        if (plot.type == "All" & !same.window) 
            dev.new()
        if (is.null(user.xlab) & !is.null(x.labels[[3]])) 
            xlab <- x.labels[[3]]
        else xlab <- user.xlab
        if (is.null(user.ylab)) {
            if (!is.null(y.labels[[3]])) 
                ylab <- y.labels[[3]]
            else {
                ylab <- paste("Quantiles of the", score.type, 
                  "Scores")
            }
        }
        else ylab <- user.ylab
        if (score.type == "Normal") {
            qqPlot(gof.obj$group.scores, distribution = "norm", 
                param.list = list(mean = 0, sd = 1), digits = digits, 
                add.line = add.line, qq.line.type = "0-1", ..., 
                cex.axis = cex.axis, cex.lab = cex.lab, xlab = user.xlab, 
                ylab = ylab, main = "", xlim = user.xlim, ylim = user.ylim)
        }
        else {
            qqPlot(gof.obj$group.scores, distribution = "chisq", 
                param.list = list(df = 2), digits = digits, add.line = add.line, 
                qq.line.type = "0-1", ..., cex.axis = cex.axis, 
                cex.lab = cex.lab, xlab = user.xlab, ylab = ylab, 
                main = "", xlim = user.xlim, ylim = user.ylim)
        }
        if (is.null(user.main)) {
            if (!is.null(captions[[3]])) 
                main <- captions[[3]]
            else {
                main <- paste("Quantile-Quantile Plot of ", score.type, 
                  " Scores for\n", paste(data.name, "Grouped by", 
                    grouping.variable), sep = "")
            }
        }
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
    }
    if (is.element(plot.type, c("All", "Scores Tukey M-D Q-Q Plot"))) {
        method <- gof.obj$method
        n.char <- nchar(method)
        score.type <- ifelse(substr(method, n.char - 13, n.char - 
            1) == "Normal Scores", "Normal", "Chi-square")
        if (plot.type == "All" & !same.window) 
            dev.new()
        if (is.null(user.xlab) & !is.null(x.labels[[4]])) 
            xlab <- x.labels[[4]]
        else xlab <- user.xlab
        if (is.null(user.ylab) & !is.null(y.labels[[4]])) 
            ylab <- y.labels[[4]]
        else ylab <- user.ylab
        if (score.type == "Normal") {
            qqPlot(gof.obj$group.scores, distribution = "norm", 
                param.list = list(mean = 0, sd = 1), plot.type = "Tukey Mean-Difference Q-Q", 
                add.line = add.line, digits = digits, ..., cex.axis = cex.axis, 
                cex.lab = cex.lab, main = "", xlab = user.xlab, 
                ylab = user.ylab, xlim = user.xlim, ylim = user.ylim)
        }
        else {
            qqPlot(gof.obj$group.scores, distribution = "chisq", 
                param.list = list(df = 2), plot.type = "Tukey Mean-Difference Q-Q", 
                add.line = add.line, digits = digits, ..., cex.axis = cex.axis, 
                cex.lab = cex.lab, main = "", xlab = user.xlab, 
                ylab = user.ylab, xlim = user.xlim, ylim = user.ylim)
        }
        if (is.null(user.main)) {
            if (!is.null(captions[[2]])) 
                main <- captions[[2]]
            else main <- paste("Tukey Mean-Difference Q-Q Plot of ", 
                score.type, " Scores for\n", paste(data.name, 
                  "Grouped by", grouping.variable), sep = "")
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
        if (!is.null(gof.obj$n.param.est) && gof.obj$n.param.est > 
            0) {
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
        text(0, 1 - mf * char.ht, "Grouped With:", adj = 0)
        text(1, 1 - mf * char.ht, grouping.variable, adj = 1)
        mf <- mf + 2
        if (!is.null(gof.obj$subset.expression)) {
            text(0, 1 - mf * char.ht, "Subset With:", adj = 0)
            text(1, 1 - mf * char.ht, gof.obj$subset.expression, 
                adj = 1)
            mf <- mf + 2
        }
        if (!is.null(gof.obj$parent.of.data)) {
            text(0, 1 - mf * char.ht, "Data Source:", adj = 0)
            text(1, 1 - mf * char.ht, gof.obj$parent.of.data, 
                adj = 1)
            mf <- mf + 2
        }
        text(0, 1 - mf * char.ht, "Number of Groups:", adj = 0)
        text(1, 1 - mf * char.ht, gof.obj$n.groups, adj = 1)
        mf <- mf + 2
        text(0, 1 - mf * char.ht, "Test Statistic:", adj = 0)
        text(1, 1 - mf * char.ht, paste(format(names(gof.obj$statistic), 
            justify = "left"), " = ", format(format(gof.obj$statistic, 
            digits = digits, nsmall = 0)), sep = "", collapse = ""), 
            adj = 1)
        mf <- mf + 2 + length(gof.obj$statistic) - 1
        if (!is.null(gof.obj$parameters)) {
            text(0, 1 - mf * char.ht, "Test Statistic Parmeter:", 
                adj = 0)
            text(1, 1 - mf * char.ht, paste(format(names(gof.obj$parameters), 
                justify = "left"), format(gof.obj$parameters, 
                digits = digits, nsmall = 0), sep = " = "), adj = 1)
            mf <- mf + 2
        }
        if (individual.p.values) {
            index <- 1:gof.obj$n.groups
            text(0, 1 - mf * char.ht, "Individual P-values:", 
                adj = c(0, 1))
            text(1, 1 - mf * char.ht, paste(format(names(gof.obj$p.value[index]), 
                justify = "left"), " = ", format(format(gof.obj$p.value[index], 
                digits = digits, nsmall = 0)), "\n", sep = "", 
                collapse = ""), adj = c(1, 1))
            mf <- mf + 2 + gof.obj$n.groups
        }
        index <- 1 + gof.obj$n.groups
        text(0, 1 - mf * char.ht, "P-value:", adj = 0)
        text(1, 1 - mf * char.ht, format(gof.obj$p.value[index], 
            digits = digits, nsmall = 0), adj = 1)
        par(c(o.mar, o.par))
    }
    if (plot.type == "Summary" & add.om.title) {
        if (is.null(om.title)) {
            string1 <- paste("Results of", gof.obj$method, "for", 
                gof.obj$distribution, "Distribution for")
            string2 <- paste(data.name.string, "Grouped with ", 
                grouping.variable)
            om.title <- paste(string1, string2, sep = "\n")
        }
        mtext(om.title, side = 3, line = om.line, outer = TRUE, 
            cex = om.cex.main, font = om.font)
    }
    invisible(gof.obj)
}
