plotTTestLnormAltDesign <-
function (x.var = "n", y.var = "power", range.x.var = NULL, n.or.n1 = 25, 
    n2 = n.or.n1, ratio.of.means = switch(alternative, greater = 2, 
        less = 0.5, two.sided = ifelse(two.sided.direction == 
            "greater", 2, 0.5)), cv = 1, alpha = 0.05, power = 0.95, 
    sample.type = ifelse(!missing(n2), "two.sample", "one.sample"), 
    alternative = "two.sided", two.sided.direction = "greater", 
    approx = FALSE, round.up = FALSE, n.max = 5000, tol = 1e-07, 
    maxiter = 1000, plot.it = TRUE, add = FALSE, n.points = 50, 
    plot.col = "black", plot.lwd = 3 * par("cex"), plot.lty = 1, 
    digits = .Options$digits, cex.main = par("cex"), ..., main = NULL, 
    xlab = NULL, ylab = NULL, type = "l") 
{
    x.var <- match.arg(x.var, c("n", "ratio.of.means", "cv", 
        "power", "alpha"))
    y.var <- match.arg(y.var, c("power", "ratio.of.means", "n"))
    sample.type <- match.arg(sample.type, c("one.sample", "two.sample"))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    two.sided.direction <- match.arg(two.sided.direction, c("greater", 
        "less"))
    if (x.var == y.var) 
        stop("'x.var' and 'y.var' cannot denote the same quantity")
    if (missing(range.x.var)) {
        range.x.var <- switch(x.var, n = c(2, 50), ratio.of.means = switch(alternative, 
            greater = c(1, 2), less = c(0.5, 1), two.sided = if (two.sided.direction == 
                "greater") c(1, 2) else c(0.5, 1)), cv = c(0.5, 
            2), power = c(alpha + .Machine$double.eps, 0.95), 
            alpha = c(0.01, 0.2))
    }
    else {
        if (is.null(range.x.var) || !all(is.finite(range.x.var)) || 
            !is.vector(range.x.var, mode = "numeric") || length(range.x.var) != 
            2) 
            stop(paste("'range.x.var' must be a numeric vector of length 2", 
                "with no missing (NA), infinite(-Inf, Inf), or undefined(NaN) values"))
    }
    min.x <- range.x.var[1]
    max.x <- range.x.var[2]
    if (min.x >= max.x) 
        stop("The second element of 'range.x.var' must be larger than the first")
    if (!is.vector(n.points, mode = "numeric") || length(n.points) != 
        1 || n.points != trunc(n.points) || n.points < 2) 
        stop("'n.points' must be an integer larger than 1")
    if (x.var != "alpha") {
        if (is.null(alpha) || !is.finite(alpha) || !is.vector(alpha, 
            mode = "numeric") || length(alpha) != 1 || alpha <= 
            0 || alpha >= 1) {
            stop("'alpha' must be a scalar between 0 and 1")
        }
    }
    switch(x.var, n = {
        if (!all(range.x.var == trunc(range.x.var)) || min.x < 
            2) stop(paste("When x.var='n', 'range.x.var' must be", 
            "a vector containing two integers, with the first", 
            "element greater than 1, and the second", "element larger than the first"))
    }, ratio.of.means = {
        if (min.x <= 0) stop(paste("When x.var='ratio.of.means', all values of", 
            "'ratio.of.means' must be positive,", "so the first element of 'range.x.var' must be positive."))
    }, cv = {
        if (min.x <= 0) stop(paste("When x.var='cv', all values of", 
            "'cv' must be positive,", "so the first element of 'range.x.var' must be positive."))
    }, power = {
        if (min.x < alpha || max.x >= 1) stop(paste("When x.var='power', 'range.x.var' must be", 
            "a vector of two elements, with the first element", 
            "greater than or equal to 'alpha', the second element", 
            "less than 1, and the second", "element larger than the first."))
    }, alpha = {
        if (min.x <= 0 || max.x >= 1) stop(paste("When x.var='alpha', 'range.x.var' must be", 
            "a vector of two elements, with the first element", 
            "greater than 0, the second element less than 1, and the second", 
            "element larger than the first."))
    })
    n2.constrained <- !missing(n2) && !is.null(n2)
    if (x.var != "n" && y.var != "n") {
        if (is.null(n.or.n1) || !is.finite(n.or.n1) || !is.vector(n.or.n1, 
            mode = "numeric") || length(n.or.n1) != 1 || n.or.n1 < 
            2 || n.or.n1 != trunc(n.or.n1)) 
            stop("'n.or.n1 must be an integer greater than 1")
        if (sample.type == "two.sample" && (is.null(n2) || !is.finite(n2) || 
            !is.vector(n2, mode = "numeric") || length(n2) != 
            1 || n2 < 2 || n2 != trunc(n2))) 
            stop("n2 must be an integer greater than 1")
    }
    else if (sample.type == "two.sample" && n2.constrained && 
        (is.null(n2) || !is.finite(n2) || !is.vector(n2, mode = "numeric") || 
            length(n2) != 1 || n2 < 2 || n2 != trunc(n2))) 
        stop("n2 must be an integer greater than 1")
    if (x.var != "ratio.of.means" && y.var != "ratio.of.means" && 
        (is.null(ratio.of.means) || !is.vector(ratio.of.means, 
            mode = "numeric") || length(ratio.of.means) != 1 || 
            ratio.of.means <= 0)) 
        stop("'ratio.of.means' must be a positive numeric scalar")
    if (x.var != "cv" && (is.null(cv) || !is.vector(cv, mode = "numeric") || 
        length(cv) != 1 || cv <= 0)) 
        stop("'cv' must be a positive numeric scalar")
    if (x.var != "power" && y.var != "power" && (is.null(power) || 
        !is.vector(power, mode = "numeric") || length(power) != 
        1 || power <= 0 || power >= 1)) 
        stop("'power' must be a numeric scalar greater than 0 and less than 1")
    if (!is.vector(n.max, mode = "numeric") || length(n.max) != 
        1 || !is.finite(n.max) || n.max != trunc(n.max) || n.max < 
        2) 
        stop("'n.max' must be a positive integer greater than 1")
    if (!is.vector(maxiter, mode = "numeric") || length(maxiter) != 
        1 || !is.finite(maxiter) || maxiter != trunc(maxiter) || 
        maxiter < 2) 
        stop("'maxiter' must be a positive integer greater than 1")
    type.string <- ifelse(sample.type == "one.sample", "One-Sample", 
        "Two-Sample")
    ratio.string <- ifelse(sample.type == "one.sample", "theta / theta0", 
        "theta1 / theta2")
    alt.string <- switch(alternative, two.sided = "(Two-Sided Alternative)", 
        less = "(Lower One-Sided Alternative)", greater = "(Upper One-Sided Alternative)")
    if (x.var != "n" && y.var != "n") {
        if (sample.type == "two.sample") 
            n.string <- paste("n1 = ", n.or.n1, ", n2 = ", n2, 
                ",", sep = "")
        else n.string <- paste("n =", n.or.n1)
    }
    if (plot.it) 
        gen.gp.list <- checkGraphicsPars(...)$gen.gp.list
    combo <- paste(sort(c(x.var, y.var)), collapse = " & ")
    switch(combo, `alpha & ratio.of.means` = {
        x <- seq(min.x, max.x, length = n.points)
        y <- tTestLnormAltRatioOfMeans(n.or.n1 = n.or.n1, n2 = n2, 
            cv = cv, alpha = x, power = power, sample.type = sample.type, 
            alternative = alternative, two.sided.direction = two.sided.direction, 
            approx = approx, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- ratio.string
        if (is.null(xlab)) xlab <- "Alpha"
        line1 <- paste("Ratio of Means vs. Alpha for ", type.string, 
            " t-Test ", alt.string, sep = "")
        line2 <- "for Lognormal Data with"
        line3 <- paste(n.string, ", CV = ", format(cv, digits = digits), 
            ", and Power = ", format(power, digits = digits), 
            sep = "")
    }, `alpha & n` = {
        if (alternative == "greater" && ratio.of.means <= 1) stop(paste("When alternative='greater',", 
            "'ratio.of.means' must be bigger than 1."))
        if (alternative == "less" && (ratio.of.means <= 0 || 
            ratio.of.means >= 1)) stop(paste("When alternative='less',", 
            "'ratio.of.means' must be between 0 and 1."))
        x <- seq(min.x, max.x, length = n.points)
        if (sample.type == "two.sample") {
            if (n2.constrained) {
                if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                  n2)
                y <- tTestLnormAltN(ratio.of.means = ratio.of.means, 
                  cv = cv, alpha = x, power = power, sample.type = "two.sample", 
                  alternative = alternative, approx = approx, 
                  n2 = n2, round.up = round.up, n.max = n.max, 
                  tol = tol, maxiter = maxiter)$n1
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                y <- tTestLnormAltN(ratio.of.means = ratio.of.means, 
                  cv = cv, alpha = x, power = power, sample.type = "two.sample", 
                  alternative = alternative, approx = approx, 
                  round.up = round.up, n.max = n.max, tol = tol, 
                  maxiter = maxiter)
            }
        } else {
            if (is.null(ylab)) ylab <- "Sample Size (n)"
            y <- tTestLnormAltN(ratio.of.means = ratio.of.means, 
                cv = cv, alpha = x, power = power, sample.type = "one.sample", 
                alternative = alternative, approx = approx, round.up = round.up, 
                n.max = n.max, tol = tol, maxiter = maxiter)
        }
        if (is.null(xlab)) xlab <- "Alpha"
        line1 <- paste("Sample Size vs. Alpha for ", type.string, 
            " t-Test ", alt.string, sep = "")
        line2 <- "for Lognormal Data with"
        line3 <- paste("Ratio of Means = ", format(ratio.of.means, 
            digits = digits), ", CV = ", format(cv, digits = digits), 
            ", and Power = ", format(power, digits = digits), 
            sep = "")
    }, `alpha & power` = {
        x <- seq(min.x, max.x, length = n.points)
        y <- tTestLnormAltPower(n.or.n1 = n.or.n1, n2 = n2, ratio.of.means = ratio.of.means, 
            cv = cv, alpha = x, sample.type = sample.type, alternative = alternative, 
            approx = approx)
        if (is.null(xlab)) xlab <- "Alpha"
        if (is.null(ylab)) ylab <- "Power"
        line1 <- paste("Power vs. Alpha for ", type.string, " t-Test ", 
            alt.string, sep = "")
        line2 <- "for Lognormal Data with"
        line3 <- paste(n.string, ", Ratio of Means = ", format(ratio.of.means, 
            digits = digits), ", and CV = ", format(cv, digits = digits), 
            sep = "")
    }, `cv & ratio.of.means` = {
        x <- seq(min.x, max.x, length = n.points)
        y <- tTestLnormAltRatioOfMeans(n.or.n1 = n.or.n1, n2 = n2, 
            cv = x, alpha = alpha, power = power, sample.type = sample.type, 
            alternative = alternative, two.sided.direction = two.sided.direction, 
            approx = approx, tol = tol, maxiter = maxiter)
        if (is.null(ylab)) ylab <- ratio.string
        if (is.null(xlab)) xlab <- "CV"
        line1 <- paste("Ratio of Means vs. CV for ", type.string, 
            " t-Test ", alt.string, sep = "")
        line2 <- "for Lognormal Data with"
        line3 <- paste(n.string, ",  Power = ", format(power, 
            digits = digits), ", and Alpha = ", format(alpha, 
            digits = digits), sep = "")
    }, `cv & n` = {
        if (alternative == "greater" && ratio.of.means <= 1) stop(paste("When alternative='greater',", 
            "'ratio.of.means' must be bigger than 1."))
        if (alternative == "less" && (ratio.of.means <= 0 || 
            ratio.of.means >= 1)) stop(paste("When alternative='less',", 
            "'ratio.of.means' must be between 0 and 1."))
        x <- seq(min.x, max.x, length = n.points)
        if (sample.type == "two.sample") {
            if (n2.constrained) {
                if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                  n2)
                y <- tTestLnormAltN(ratio.of.means = ratio.of.means, 
                  cv = x, alpha = alpha, power = power, sample.type = "two.sample", 
                  alternative = alternative, approx = approx, 
                  n2 = n2, round.up = round.up, n.max = n.max, 
                  tol = tol, maxiter = maxiter)$n1
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                y <- tTestLnormAltN(ratio.of.means = ratio.of.means, 
                  cv = x, alpha = alpha, power = power, sample.type = "two.sample", 
                  alternative = alternative, approx = approx, 
                  round.up = round.up, n.max = n.max, tol = tol, 
                  maxiter = maxiter)
            }
        } else {
            if (is.null(ylab)) ylab <- "Sample Size (n)"
            y <- tTestLnormAltN(ratio.of.means = ratio.of.means, 
                cv = x, alpha = alpha, power = power, sample.type = "one.sample", 
                alternative = alternative, approx = approx, round.up = round.up, 
                n.max = n.max, tol = tol, maxiter = maxiter)
        }
        if (is.null(xlab)) xlab <- "CV"
        line1 <- paste("Sample Size vs. CV for ", type.string, 
            " t-Test ", alt.string, sep = "")
        line2 <- "for Lognormal Data with"
        line3 <- paste("Ratio of Means = ", format(ratio.of.means, 
            digits = digits), ", Power = ", format(power, digits = digits), 
            ", and Alpha = ", format(alpha, digits = digits), 
            sep = "")
    }, `cv & power` = {
        x <- seq(min.x, max.x, length = n.points)
        y <- tTestLnormAltPower(n.or.n1 = n.or.n1, n2 = n2, ratio.of.means = ratio.of.means, 
            cv = x, alpha = alpha, sample.type = sample.type, 
            alternative = alternative, approx = approx)
        if (is.null(xlab)) xlab <- "CV"
        if (is.null(ylab)) ylab <- "Power"
        line1 <- paste("Power vs. CV for ", type.string, " t-Test ", 
            alt.string, sep = "")
        line2 <- "for Lognormal Data with"
        line3 <- paste(n.string, ", Ratio of Means = ", format(ratio.of.means, 
            digits = digits), ", and Alpha = ", format(alpha, 
            digits = digits), sep = "")
    }, `n & ratio.of.means` = {
        if (x.var == "ratio.of.means") {
            if (alternative == "greater" && min.x <= 1) stop(paste("When alternative='greater',", 
                "and x.var='ratio.of.means' and y.var='n', all values of", 
                "'ratio.of.means' must be bigger than 1,", "so the first element of 'range.x.var' must be bigger than 1."))
            if (alternative == "less" && (min.x <= 0 || max.x >= 
                1)) stop(paste("When alternative='less',", "and x.var='ratio.of.means' and y.var='n', all values of", 
                "'ratio.of.means' must be between 0 and 1,", 
                "so the first element of 'range.x.var' must be greater than 0", 
                "and the second element of 'range.x.var' must be less than 1."))
            x <- seq(min.x, max.x, length = n.points)
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                  y <- tTestLnormAltN(ratio.of.means = x, cv = cv, 
                    alpha = alpha, power = power, sample.type = "two.sample", 
                    alternative = alternative, approx = approx, 
                    n2 = n2, round.up = round.up, n.max = n.max, 
                    tol = tol, maxiter = maxiter)$n1
                } else {
                  if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                  y <- tTestLnormAltN(ratio.of.means = x, cv = cv, 
                    alpha = alpha, power = power, sample.type = "two.sample", 
                    alternative = alternative, approx = approx, 
                    round.up = round.up, n.max = n.max, tol = tol, 
                    maxiter = maxiter)
                }
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n)"
                y <- tTestLnormAltN(ratio.of.means = x, cv = cv, 
                  alpha = alpha, power = power, sample.type = "one.sample", 
                  alternative = alternative, approx = approx, 
                  round.up = round.up, n.max = n.max, tol = tol, 
                  maxiter = maxiter)
            }
            if (is.null(xlab)) xlab <- ratio.string
            line1 <- paste("Sample Size vs. Ratio of Means for ", 
                type.string, " t-Test ", alt.string, sep = "")
            line2 <- "for Lognormal Data with"
            line3 <- paste("CV = ", format(cv, digits = digits), 
                ", Power = ", format(power, digits = digits), 
                ", and Alpha = ", format(alpha, digits = digits), 
                sep = "")
        } else {
            x <- seq(min.x, max.x, length = n.points)
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(xlab)) xlab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                } else {
                  n2 <- x
                  if (is.null(xlab)) xlab <- "Sample Size (n1 and n2)"
                }
            } else if (is.null(xlab)) xlab <- "Sample Size (n)"
            y <- tTestLnormAltRatioOfMeans(n.or.n1 = x, n2 = n2, 
                cv = cv, alpha = alpha, power = power, sample.type = sample.type, 
                alternative = alternative, two.sided.direction = two.sided.direction, 
                approx = approx, tol = tol, maxiter = maxiter)
            if (is.null(ylab)) ylab <- ratio.string
            line1 <- paste("Ratio of Means vs. Sample Size for ", 
                type.string, " t-Test ", alt.string, sep = "")
            line2 <- "for Lognormal Data with"
            line3 <- paste("CV = ", format(cv, digits = digits), 
                ", Power = ", format(power, digits = digits), 
                ", and Alpha = ", format(alpha, digits = digits), 
                sep = "")
        }
    }, `power & ratio.of.means` = {
        if (x.var == "ratio.of.means") {
            x <- seq(min.x, max.x, length = n.points)
            y <- tTestLnormAltPower(n.or.n1 = n.or.n1, n2 = n2, 
                ratio.of.means = x, cv = cv, alpha = alpha, sample.type = sample.type, 
                alternative = alternative, approx = approx)
            if (is.null(xlab)) xlab <- ratio.string
            if (is.null(ylab)) ylab <- "Power"
            line1 <- paste("Power vs. Ratio of Means for ", type.string, 
                " t-Test ", alt.string, sep = "")
            line2 <- "for Lognormal Data with"
            line3 <- paste("CV = ", format(cv, digits = digits), 
                ", ", n.string, ", and Alpha = ", format(alpha, 
                  digits = digits), sep = "")
        } else {
            x <- seq(min.x, max.x, length = n.points)
            y <- tTestLnormAltRatioOfMeans(n.or.n1 = n.or.n1, 
                n2 = n2, cv = cv, alpha = alpha, power = x, sample.type = sample.type, 
                alternative = alternative, two.sided.direction = two.sided.direction, 
                approx = approx, tol = tol, maxiter = maxiter)
            if (is.null(ylab)) ylab <- ratio.string
            if (is.null(xlab)) xlab <- "Power"
            line1 <- paste("Ratio of Means vs. Power for ", type.string, 
                " t-Test ", alt.string, sep = "")
            line2 <- "for Lognormal Data with"
            line3 <- paste("CV = ", format(cv, digits = digits), 
                ", ", n.string, ", and Alpha = ", format(alpha, 
                  digits = digits), sep = "")
        }
    }, `n & power` = {
        if (x.var == "n") {
            x <- seq(min.x, max.x, length = n.points)
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(xlab)) xlab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                } else {
                  n2 <- x
                  if (is.null(xlab)) xlab <- "Sample Size (n1 and n2)"
                }
            } else if (is.null(xlab)) xlab <- "Sample Size (n)"
            y <- tTestLnormAltPower(n.or.n1 = x, n2 = n2, ratio.of.means = ratio.of.means, 
                cv = cv, alpha = alpha, sample.type = sample.type, 
                alternative = alternative, approx = approx)
            if (is.null(ylab)) ylab <- "Power"
            line1 <- paste("Power vs. Sample Size for ", type.string, 
                " t-Test ", alt.string, sep = "")
            line2 <- "for Lognormal Data with"
            line3 <- paste("Ratio of Means = ", format(ratio.of.means, 
                digits = digits), ", CV = ", format(cv, digits = digits), 
                ", and Alpha = ", format(alpha, digits = digits), 
                sep = "")
        } else {
            if (alternative == "greater" && ratio.of.means <= 
                1) stop(paste("When alternative='greater',", 
                "'ratio.of.means' must be bigger than 1."))
            if (alternative == "less" && (ratio.of.means <= 0 || 
                ratio.of.means >= 1)) stop(paste("When alternative='less',", 
                "'ratio.of.means' must be between 0 and 1."))
            x <- seq(min.x, max.x, length = n.points)
            if (sample.type == "two.sample") {
                if (n2.constrained) {
                  if (is.null(ylab)) ylab <- paste("Sample Size (n1), With n2 =", 
                    n2)
                  y <- tTestLnormAltN(ratio.of.means = ratio.of.means, 
                    cv = cv, alpha = alpha, power = x, sample.type = "two.sample", 
                    alternative = alternative, approx = approx, 
                    n2 = n2, round.up = round.up, n.max = n.max, 
                    tol = tol, maxiter = maxiter)$n1
                } else {
                  if (is.null(ylab)) ylab <- "Sample Size (n1 and n2)"
                  y <- tTestLnormAltN(ratio.of.means = ratio.of.means, 
                    cv = cv, alpha = alpha, power = x, sample.type = "two.sample", 
                    alternative = alternative, approx = approx, 
                    round.up = round.up, n.max = n.max, tol = tol, 
                    maxiter = maxiter)
                }
            } else {
                if (is.null(ylab)) ylab <- "Sample Size (n)"
                y <- tTestLnormAltN(ratio.of.means = ratio.of.means, 
                  cv = cv, alpha = alpha, power = x, sample.type = "one.sample", 
                  alternative = alternative, approx = approx, 
                  round.up = round.up, n.max = n.max, tol = tol, 
                  maxiter = maxiter)
            }
            if (is.null(xlab)) xlab <- "Power"
            line1 <- paste("Sample Size vs. Power for ", type.string, 
                " t-Test ", alt.string, sep = "")
            line2 <- "for Lognormal Data with"
            line3 <- paste("Ratio of Means = ", format(ratio.of.means, 
                digits = digits), ", CV = ", format(cv, digits = digits), 
                ", and Alpha =", format(alpha, digits = digits), 
                sep = "")
        }
    })
    if (plot.it) {
        if (!add) {
            plot(x, y, type = "n", main = "", sub = "", ..., 
                xlab = xlab, ylab = ylab)
            if (is.null(main)) {
                mtext(text = line1, side = 3, line = 3, cex = cex.main)
                mtext(text = line2, side = 3, line = 2, cex = cex.main)
                mtext(text = line3, side = 3, line = 1, cex = cex.main)
            }
            else {
                arg.list <- c(list(main = main), gen.gp.list, 
                  list(cex = cex.main))
                do.call("title", arg.list)
            }
            arg.list <- c(list(x = x, y = y), gen.gp.list, list(type = type, 
                col = plot.col, lwd = plot.lwd, lty = plot.lty))
            do.call("lines", arg.list)
        }
        else {
            arg.list <- c(list(x = x, y = y), gen.gp.list, list(type = type, 
                col = plot.col, lwd = plot.lwd, lty = plot.lty))
            do.call("lines", arg.list)
        }
    }
    ret.list <- list(x, y)
    names(ret.list) <- c(x.var, y.var)
    invisible(ret.list)
}
