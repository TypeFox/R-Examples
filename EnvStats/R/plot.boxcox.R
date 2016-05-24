plot.boxcox <-
function (x, plot.type = "Objective vs. lambda", same.window = TRUE, 
    ask = same.window & plot.type != "Ojective vs. lambda", plot.pos.con = 0.375, 
    estimate.params = FALSE, equal.axes = qq.line.type == "0-1" || 
        estimate.params, add.line = TRUE, qq.line.type = "least squares", 
    duplicate.points.method = "standard", points.col = 1, line.col = 1, 
    line.lwd = par("cex"), line.lty = 1, digits = .Options$digits, 
    cex.main = 1.4 * par("cex"), cex.sub = par("cex"), main = NULL, 
    sub = NULL, xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, 
    ...) 
{
    boxcox.obj <- x
    plot.type <- match.arg(plot.type, c("All", "Objective vs. lambda", 
        "Q-Q Plots", "Tukey M-D Q-Q Plots"))
    qq.line.type <- match.arg(qq.line.type, c("least squares", 
        "0-1", "robust"))
    duplicate.points.method <- match.arg(duplicate.points.method, 
        c("standard", "jitter", "number"))
    lambda <- boxcox.obj$lambda
    objective <- boxcox.obj$objective
    objective.name <- boxcox.obj$objective.name
    check.gp.list <- checkGraphicsPars(...)
    gp.arg.list <- check.gp.list$gp.arg.list
    gen.gp.list <- check.gp.list$gen.gp.list
    data.name <- boxcox.obj$data.name
    data.name.string <- data.name
    parent.of.data <- boxcox.obj$parent.of.data
    if (!is.null(parent.of.data)) 
        data.name.string <- paste(data.name, "in", parent.of.data)
    user.main <- main
    user.sub <- sub
    user.xlab <- xlab
    user.ylab <- ylab
    user.xlim <- xlim
    user.ylim <- ylim
    if (plot.type != "Objective vs. lambda" & same.window) {
        devAskNewPage(ask = ask)
    }
    if (is.element(plot.type, c("All", "Objective vs. lambda"))) {
        if (is.null(user.xlab)) 
            xlab <- expression(paste(lambda))
        if (is.null(user.ylab)) 
            ylab <- objective.name
        arg.list <- list(x = lambda, y = objective, type = "n", 
            xlab = xlab, ylab = ylab)
        arg.list <- c(arg.list, ...)
        if (!is.null(user.xlim)) 
            arg.list <- c(arg.list, list(xlim = user.xlim))
        if (!is.null(user.ylim)) 
            arg.list <- c(arg.list, list(ylim = user.ylim))
        do.call("plot", arg.list)
        arg.list <- c(list(x = lambda, y = objective, col = points.col), 
            gen.gp.list)
        do.call("points", arg.list)
        if (is.null(user.main)) {
            string1 <- "Box-Cox Transformation Results:"
            string2 <- paste(objective.name, " vs. lambda for ", 
                data.name.string, sep = "")
            main <- paste(string1, string2, sep = "\n")
        }
        else main <- user.main
        arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
            list(main = main))
        do.call("title", arg.list)
    }
    if (is.element(plot.type, c("All", "Q-Q Plots"))) {
        data <- boxcox.obj$data
        if (is.null(data)) 
            stop("'data' component missing from 'x'")
        eps <- boxcox.obj$eps
        n <- length(lambda)
        for (i in 1:n) {
            lambda.i <- lambda[i]
            if (is.null(user.ylab)) {
                if (lambda.i == 0) 
                  qlab <- paste("log(", data.name, ")")
                else qlab <- paste("[ (", data.name, "^", lambda[i], 
                  "- 1 )", "/", lambda[i], "]")
                ylab <- paste("Quantiles of", qlab)
            }
            y <- boxcoxTransform(data, lambda = lambda[i], eps = eps)
            if (!same.window) 
                dev.new()
            qqPlot(y, estimate.params = estimate.params, plot.pos.con = plot.pos.con, 
                equal.axes = equal.axes, add.line = add.line, 
                qq.line.type = qq.line.type, duplicate.points.method = duplicate.points.method, 
                points.col = points.col, line.col = line.col, 
                line.lwd = line.lwd, line.lty = line.lty, digits = digits, 
                ..., xlab = user.xlab, ylab = ylab, main = "", 
                xlim = user.xlim, ylim = user.ylim)
            if (is.null(user.main)) {
                string1 <- "Normal Q-Q Plot of Box-Cox Transformation for"
                string2 <- paste(data.name.string, " with lambda = ", 
                  format(lambda[i], digits = digits), sep = "")
                main <- paste(string1, string2, sep = "\n")
            }
            else main <- user.main
            arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
                list(main = main))
            do.call("title", arg.list)
            if (is.null(user.sub)) {
                sub <- paste(objective.name, "=", format(objective[i], 
                  digits = digits))
            }
            else sub <- user.sub
            mtext(sub, side = 1, line = 4, cex = cex.sub)
        }
    }
    if (is.element(plot.type, c("All", "Tukey M-D Q-Q Plots"))) {
        data <- boxcox.obj$data
        if (is.null(data)) 
            stop("'data' component missing from 'x'")
        eps <- boxcox.obj$eps
        n <- length(lambda)
        for (i in 1:n) {
            y <- boxcoxTransform(data, lambda = lambda[i], eps = eps)
            if (!same.window) 
                dev.new()
            qqPlot(y, estimate.params = TRUE, plot.type = "Tukey Mean-Difference Q-Q", 
                plot.pos.con = plot.pos.con, add.line = add.line, 
                duplicate.points.method = duplicate.points.method, 
                points.col = points.col, line.col = line.col, 
                line.lwd = line.lwd, line.lty = line.lty, digits = digits, 
                ..., main = "", xlab = user.xlab, ylab = user.ylab, 
                xlim = user.xlim, ylim = user.ylim)
            if (is.null(user.main)) {
                string1 <- "Mean-Difference Q-Q Plot of Box-Cox Transformation for"
                string2 <- paste(data.name.string, " with lambda = ", 
                  format(lambda[i], digits = digits), sep = "")
                main <- paste(string1, string2, sep = "\n")
            }
            else main <- user.main
            arg.list <- c(list(cex.main = cex.main), gen.gp.list, 
                list(main = main))
            do.call("title", arg.list)
            if (is.null(user.sub)) {
                sub <- paste(objective.name, "=", format(objective[i], 
                  digits = digits))
            }
            else sub <- user.sub
            mtext(sub, side = 1, line = 4, cex = cex.sub)
        }
    }
    invisible(boxcox.obj)
}
