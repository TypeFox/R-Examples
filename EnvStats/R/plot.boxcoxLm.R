plot.boxcoxLm <-
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
    boxcoxLm.obj <- x
    plot.type <- match.arg(plot.type, c("All", "Objective vs. lambda", 
        "Q-Q Plots", "Tukey M-D Q-Q Plots"))
    qq.line.type <- match.arg(qq.line.type, c("least squares", 
        "0-1", "robust"))
    duplicate.points.method <- match.arg(duplicate.points.method, 
        c("standard", "jitter", "number"))
    lambda <- boxcoxLm.obj$lambda
    objective <- boxcoxLm.obj$objective
    objective.name <- boxcoxLm.obj$objective.name
    check.gp.list <- checkGraphicsPars(...)
    gp.arg.list <- check.gp.list$gp.arg.list
    gen.gp.list <- check.gp.list$gen.gp.list
    data.name <- boxcoxLm.obj$data.name
    data.name.string <- data.name
    lm.obj <- boxcoxLm.obj$lm.obj
    response.name <- deparse(formula(lm.obj)[[2]])
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
        lm.obj <- boxcoxLm.obj$lm.obj
        if (is.null(lm.obj$y) || is.null(lm.obj$qr)) 
            lm.obj <- update(lm.obj, y = TRUE, qr = TRUE)
        Y <- lm.obj$y
        eps <- boxcoxLm.obj$eps
        n <- length(lambda)
        for (i in 1:n) {
            lambda.i <- lambda[i]
            if (is.null(user.ylab)) {
                if (lambda.i == 0) 
                  qlab <- paste("log(", response.name, ")")
                else qlab <- paste("[ (", response.name, "^", 
                  lambda.i, "- 1 )", "/", lambda.i, "]")
                ylab <- paste("Quantiles of Residuals Based on", 
                  qlab)
            }
            new.Y <- boxcoxTransform(x = Y, lambda = lambda[i], 
                eps = eps)
            y <- qr.resid(lm.obj$qr, new.Y)
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
        lm.obj <- boxcoxLm.obj$lm.obj
        if (is.null(lm.obj$y) || is.null(lm.obj$qr)) 
            lm.obj <- update(lm.obj, y = TRUE, qr = TRUE)
        Y <- lm.obj$y
        eps <- boxcoxLm.obj$eps
        n <- length(lambda)
        for (i in 1:n) {
            new.Y <- boxcoxTransform(x = Y, lambda = lambda[i], 
                eps = eps)
            y <- qr.resid(lm.obj$qr, new.Y)
            if (!same.window) 
                dev.new()
            qqPlot(y, estimate.params = TRUE, plot.type = "Tukey Mean-Difference Q-Q", 
                plot.pos.con = plot.pos.con, add.line = add.line, 
                duplicate.points.method = duplicate.points.method, 
                points.col = points.col, line.col = line.col, 
                line.lwd = line.lwd, line.lty = line.lty, digits = digits, 
                ..., xlab = user.xlab, ylab = user.ylab, main = "", 
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
    invisible(boxcoxLm.obj)
}
