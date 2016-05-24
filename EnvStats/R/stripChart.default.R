stripChart.default <-
function (x, method = ifelse(paired && paired.lines, "overplot", 
    "stack"), seed = 47, jitter = 0.1 * cex, offset = 1/2, vertical = TRUE, 
    group.names, group.names.cex = cex, drop.unused.levels = TRUE, 
    add = FALSE, at = NULL, xlim = NULL, ylim = NULL, ylab = NULL, 
    xlab = NULL, dlab = "", glab = "", log = "", pch = 1, col = par("fg"), 
    cex = par("cex"), points.cex = cex, axes = TRUE, frame.plot = axes, 
    show.ci = TRUE, location.pch = 16, location.cex = cex, conf.level = 0.95, 
    min.n.for.ci = 2, ci.offset = 3/ifelse(n > 2, (n - 1)^(1/3), 
        1), ci.bar.lwd = cex, ci.bar.ends = TRUE, ci.bar.ends.size = 0.5 * 
        cex, ci.bar.gap = FALSE, n.text = "bottom", n.text.line = ifelse(n.text == 
        "bottom", 2, 0), n.text.cex = cex, location.scale.text = "top", 
    location.scale.digits = 1, nsmall = location.scale.digits, 
    location.scale.text.line = ifelse(location.scale.text == 
        "top", 0, 3.5), location.scale.text.cex = cex * 0.8 * 
        ifelse(n > 6, max(0.4, 1 - (n - 6) * 0.06), 1), p.value = FALSE, 
    p.value.digits = 3, p.value.line = 2, p.value.cex = cex, 
    group.difference.ci = p.value, group.difference.conf.level = 0.95, 
    group.difference.digits = location.scale.digits, ci.and.test = "parametric", 
    ci.arg.list = NULL, test.arg.list = NULL, alternative = "two.sided", 
    plot.diff = FALSE, diff.col = col[1], diff.pch = pch[1], 
    paired = FALSE, paired.lines = paired, paired.lty = 1:6, 
    paired.lwd = 1, paired.pch = 1:14, paired.col = NULL, diff.name = NULL, 
    diff.name.cex = group.names.cex, sep.line = TRUE, sep.lty = 2, 
    sep.lwd = cex, sep.col = "gray", diff.lim = NULL, diff.at = NULL, 
    diff.axis.label = NULL, plot.diff.mar = c(5, 4, 4, 4) + 0.1, 
    ...) 
{
    method <- pmatch(method, c("overplot", "jitter", "stack"))[1L]
    if (is.na(method) || method == 0L) 
        stop("invalid plotting method")
    if (method == 2L) 
        set.seed(seed)
    n.text <- match.arg(n.text, c("bottom", "top", "none"))
    if (!vertical & n.text != "none") 
        stop("You must set vertical=T when n.text != \"none\"")
    location.scale.text <- match.arg(location.scale.text, c("none", 
        "top", "bottom"))
    if (!vertical & location.scale.text != "none") 
        stop("You must set vertical=T when location.scale.text != \"none\"")
    if (show.ci && !(length(conf.level) == 1 & is.numeric(conf.level) & 
        conf.level > 0 & conf.level < 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if (min.n.for.ci != round(min.n.for.ci) || min.n.for.ci < 
        2) 
        stop("'min.n.for.ci' must be an integer greater than 1")
    if (is.list(x)) {
        if (!all(sapply(x, is.numeric))) 
            stop("When 'x' is a data frame or list, all components must be numeric")
        groups <- x
    }
    else if (is.matrix(x)) {
        if (!is.numeric(x)) 
            stop("When 'x' is a matrix it must be numeric")
        groups <- data.frame(x)
    }
    else if (is.vector(x)) {
        if (!is.numeric(x)) 
            stop("When 'x' is a vector it must be numeric")
        groups <- list(x)
    }
    else stop("'x' must be a list, data frame, matrix, or vector")
    if (0L == (n <- length(groups))) 
        stop("invalid first argument")
    if (n == 2 && paired) {
        g1 <- groups[[1]]
        g2 <- groups[[2]]
        if (length(g1) != length(g2)) {
            stop(paste("When there are two groups and paired=TRUE,", 
                "the two groups must have the same number of observations"))
        }
        na.index <- is.na(g1) | !is.finite(g1) | is.na(g2) | 
            !is.finite(g2)
        if (any(na.index)) {
            if (sum(na.index) == length(na.index)) {
                stop("No paired observations where both values are non-missing and finite")
            }
            else {
                groups[[1]] <- g1[!na.index]
                groups[[2]] <- g2[!na.index]
            }
        }
    }
    else {
        if (any(sapply(groups, function(x) any(is.na(x)) | any(!is.finite(x))))) 
            groups <- lapply(groups, function(x) x[!is.na(x) & 
                is.finite(x)])
        groups.N <- sapply(groups, length)
        if (all(groups.N == 0)) 
            stop("All observations are missing or not finite.")
        if (drop.unused.levels) {
            groups.N <- sapply(groups, length)
            groups <- groups[groups.N > 0]
            n <- length(groups)
        }
    }
    if (n == 2L && paired && paired.lines && method == 3L) {
        stop(paste("When there are two groups and paired=TRUE and paired.lines=TRUE,", 
            "you must set method='overplot' or method='jitter'."))
    }
    if (!missing(group.names)) {
        group.names <- as.character(group.names)
        if (length(group.names) != n) 
            stop(paste("'group.names' must have the same number of elements as the number of groups.", 
                "Note that drop.unused.levels =", drop.unused.levels))
        attr(groups, "names") <- group.names
    }
    else if (is.null(attr(groups, "names"))) 
        attr(groups, "names") <- 1L:n
    group.names <- names(groups)
    if (is.null(at)) {
        if (n == 2L && plot.diff) 
            at <- 1L:3L
        else at <- 1L:n
    }
    else {
        if (n == 2L && plot.diff) {
            if (length(at) != 3) 
                stop(paste("'at' must have length equal to 3 when there are 2 groups and", 
                  "plot.diff=TRUE"))
        }
        else {
            if (length(at) != n) 
                stop(paste("'at' must have length equal to the number of groups, which is", 
                  n))
        }
    }
    if (show.ci) {
        length.ci.offset <- length(ci.offset)
        if (n == 2L && plot.diff) {
            if (!is.vector(ci.offset, mode = "numeric") || !all(is.finite(ci.offset)) || 
                any(ci.offset < .Machine$double.eps) || !(length.ci.offset %in% 
                c(1, 3))) 
                stop(paste("When there are 2 groups, plot.diff=TRUE, and show.ci=TRUE,", 
                  "'ci.offset' must be a positive scalar, or else a vector of", 
                  "positive numbers with length equal to 3"))
            if (length.ci.offset == 1) 
                ci.offset <- rep(ci.offset, 3)
        }
        else {
            if (!is.vector(ci.offset, mode = "numeric") || !all(is.finite(ci.offset)) || 
                any(ci.offset < .Machine$double.eps) || !(length.ci.offset %in% 
                c(1, n))) 
                stop(paste("'ci.offset' must be a positive scalar, or else a vector of", 
                  "positive numbers with length equal to the number of groups"))
            if (length.ci.offset == 1) 
                ci.offset <- rep(ci.offset, n)
        }
    }
    if (is.null(dlab)) 
        dlab <- deparse(substitute(x))
    if (!add) {
        if (n == 2 && plot.diff) 
            par(mar = plot.diff.mar)
        dlim <- c(NA, NA)
        dlim <- range(sapply(groups, function(x) {
            range(x[is.finite(x)], na.rm = TRUE)
        }), na.rm = TRUE)
        glim <- c(0.625, n)
        if (n == 2 && plot.diff) 
            glim[2] <- 3
        if (method == 1L & show.ci) 
            glim <- glim + if (n == 1L) 
                c(-1, 1)
            else c(0, 0.5)
        if (method == 2L) {
            glim <- glim + jitter * if (n == 1L) 
                c(-5, 5)
            else c(-2, 2)
        }
        else if (method == 3L) {
            glim <- glim + if (n == 1L) 
                c(-1, 1)
            else c(0, 0.5)
        }
        if (is.null(xlim)) 
            xlim <- if (vertical) 
                glim
            else dlim
        if (is.null(ylim)) 
            ylim <- if (vertical) 
                dlim
            else glim
        plot(xlim, ylim, type = "n", ann = FALSE, axes = FALSE, 
            log = log, ...)
        if (frame.plot) 
            box(...)
        if (axes && n == 2 && plot.diff && is.null(diff.name)) 
            diff.name <- paste(rev(group.names), collapse = "-")
        if (vertical) {
            if (axes) {
                if (n == 2L && plot.diff) {
                  axis(1, at = at[1:2], labels = group.names, 
                    cex.axis = group.names.cex, ...)
                  axis(1, at = at[3], labels = diff.name, cex.axis = diff.name.cex, 
                    ...)
                }
                else {
                  if (n > 1L) 
                    axis(1, at = at, labels = group.names, cex.axis = group.names.cex, 
                      ...)
                }
                Axis(x, side = 2, ...)
            }
            if (is.null(ylab)) 
                ylab <- dlab
            if (is.null(xlab)) 
                xlab <- glab
        }
        else {
            if (axes) {
                Axis(x, side = 1, ...)
                if (n == 2L && plot.diff) {
                  axis(2, at = at[1:2], labels = group.names, 
                    cex.axis = group.names.cex, ...)
                  axis(2, at = at[3], labels = diff.name, cex.axis = diff.name.cex, 
                    ...)
                }
                else {
                  if (n > 1L) 
                    axis(2, at = at, labels = group.names, cex.axis = group.names.cex, 
                      ...)
                }
            }
            if (is.null(xlab)) 
                xlab <- dlab
            if (is.null(ylab)) 
                ylab <- glab
        }
        title(xlab = xlab, ylab = ylab, ...)
    }
    csize <- cex * ifelse(vertical, xinch(par("cin")[1L]), yinch(par("cin")[2L]))
    n.vec <- rep.int(as.numeric(NA), n)
    location.vec <- rep.int(as.numeric(NA), n)
    scale.vec <- location.vec
    ci.mat <- matrix(as.numeric(NA), nrow = n, ncol = 3)
    ci.and.test <- match.arg(ci.and.test, c("parametric", "nonparametric"))
    if (ci.and.test == "parametric") {
        ci.fcn <- "t.test"
        location.fcn <- "mean"
        scale.fcn <- "sd"
    }
    else {
        ci.fcn <- "wilcox.test"
        location.fcn <- "median"
        scale.fcn <- "iqr"
    }
    if (ci.and.test == "nonparametric") {
        if (is.null(ci.arg.list) || all(is.na(pmatch(names(ci.arg.list), 
            "conf.int")))) {
            ci.arg.list <- c(ci.arg.list, list(conf.int = TRUE))
        }
        else {
            index <- (1:length(ci.arg.list))[!is.na(pmatch(names(ci.arg.list), 
                "conf.int"))]
            if (!unlist(ci.arg.list[index])) 
                ci.arg.list[[index]] <- TRUE
        }
    }
    dimnames(ci.mat) <- list(group.names, c("LCL", "UCL", "Conf.Level"))
    if (n == 2 && paired && paired.lines) {
        y1.at <- at[1] + (ci.offset[1] * csize)/2
        y2.at <- at[2] - (ci.offset[2] * csize)/2
        x1 <- groups[[1]]
        x2 <- groups[[2]]
        y1 <- rep.int(y1.at, length(x1))
        y2 <- rep.int(y2.at, length(x2))
        if (method == 2L) {
            y1 <- y1 + stats::runif(length(y1), -jitter, jitter)
            y2 <- y2 + stats::runif(length(y2), -jitter, jitter)
        }
        else if (method == 3L) {
            xg <- split(x1, factor(x1))
            xo <- lapply(xg, seq_along)
            x1 <- unlist(xg, use.names = FALSE)
            y1 <- rep.int(y1.at, length(x1)) + (unlist(xo, use.names = FALSE) - 
                1) * offset * csize
            xg <- split(x2, factor(x2))
            xo <- lapply(xg, seq_along)
            x2 <- unlist(xg, use.names = FALSE)
            y2 <- rep.int(y2.at, length(x2)) + (unlist(xo, use.names = FALSE) - 
                1) * offset * csize
        }
        if (is.null(paired.col)) 
            paired.col <- c("black", "red", "green3", "blue", 
                "magenta", "darkgreen", "purple", "orange", "darkolivegreen", 
                "steelblue", "darkgray")
        if (vertical) {
            points(y1, x1, col = paired.col, pch = paired.pch, 
                cex = points.cex, ...)
            points(y2, x2, col = paired.col, pch = paired.pch, 
                cex = points.cex, ...)
            segments(x0 = y1, y0 = x1, x1 = y2, y1 = x2, col = paired.col, 
                lty = paired.lty, lwd = paired.lwd, ...)
        }
        else {
            points(x1, y1, col = paired.col, pch = paired.pch, 
                cex = points.cex, ...)
            points(x2, y2, col = paired.col, pch = paired.pch, 
                cex = points.cex, ...)
            segments(x0 = x1, y0 = y1, x1 = x2, y1 = y2, col = paired.col, 
                lty = paired.lty, lwd = paired.lwd, ...)
        }
    }
    for (i in 1L:n) {
        x <- groups[[i]]
        if (!(n == 2 && paired && paired.lines)) {
            y <- rep.int(at[i], length(x))
            if (method == 2L) 
                y <- y + stats::runif(length(y), -jitter, jitter)
            else if (method == 3L) {
                xg <- split(x, factor(x))
                xo <- lapply(xg, seq_along)
                x <- unlist(xg, use.names = FALSE)
                y <- rep.int(at[i], length(x)) + (unlist(xo, 
                  use.names = FALSE) - 1) * offset * csize
            }
            if (vertical) {
                points(y, x, col = col[(i - 1L)%%length(col) + 
                  1L], pch = pch[(i - 1L)%%length(pch) + 1L], 
                  cex = points.cex, ...)
            }
            else {
                points(x, y, col = col[(i - 1L)%%length(col) + 
                  1L], pch = pch[(i - 1L)%%length(pch) + 1L], 
                  cex = points.cex, ...)
            }
        }
        x <- x[is.finite(x)]
        n.vec[i] <- sum(is.finite(x))
        if (n.text != "none") {
            side <- ifelse(n.text == "top", 3, 1)
            mtext(paste("n=", n.vec[i], sep = ""), side = side, 
                line = n.text.line, at = at[i], cex = n.text.cex)
        }
        if (n.vec[i] == 0) 
            next
        location.vec[i] <- do.call(location.fcn, args = list(x = x))
        scale.vec[i] <- do.call(scale.fcn, args = list(x = x))
        if (location.scale.text != "none") {
            side <- ifelse(location.scale.text == "top", 3, 1)
            dum <- format(round(c(location.vec[i], scale.vec[i]), 
                location.scale.digits), nsmall = nsmall)
            if (ci.and.test == "parametric") {
                string1 <- "Mean="
                string2 <- "\nSD   ="
            }
            else {
                string1 <- "Median="
                string2 <- "\nIQR     ="
            }
            mtext(paste(string1, dum[1], string2, dum[2], sep = ""), 
                side = side, line = location.scale.text.line, 
                at = at[i], cex = location.scale.text.cex)
        }
        if (n.vec[i] >= 2) {
            sd.x <- sd(x)
            if (sd.x > 0) 
                ci.vec <- do.call(ci.fcn, args = c(list(x = x, 
                  conf.level = conf.level), ci.arg.list))$conf.int
            else ci.vec <- c(NA, NA)
            ci.mat[i, c("LCL", "UCL")] <- ci.vec
            ci.mat[i, "Conf.Level"] <- conf.level
        }
        if (show.ci) {
            if (n == 2 && paired && paired.lines && i == 1) {
                if (vertical) {
                  if (sd.x > 0 & n.vec[i] >= min.n.for.ci) 
                    errorBar(x = at[i] - (ci.offset[i] * csize)/2, 
                      y = location.vec[i], lower = ci.mat[i, 
                        "LCL"], upper = ci.mat[i, "UCL"], incr = FALSE, 
                      bar.ends = ci.bar.ends, gap = ci.bar.gap, 
                      add = TRUE, horizontal = FALSE, col = col[(i - 
                        1L)%%length(col) + 1L], lwd = ci.bar.lwd, 
                      bar.ends.size = ci.bar.ends.size)
                  points(at[i] - (ci.offset[i] * csize)/2, location.vec[i], 
                    pch = location.pch, col = col[(i - 1L)%%length(col) + 
                      1L], cex = location.cex)
                }
                else {
                  if (sd.x > 0 & n.vec[i] >= min.n.for.ci) 
                    errorBar(x = location.vec[i], y = at[i] - 
                      (ci.offset[i] * csize)/2, lower = ci.mat[i, 
                      "LCL"], upper = ci.mat[i, "UCL"], incr = FALSE, 
                      bar.ends = ci.bar.ends, gap = ci.bar.gap, 
                      add = TRUE, horizontal = TRUE, col = col[(i - 
                        1L)%%length(col) + 1L], lwd = ci.bar.lwd, 
                      bar.ends.size = ci.bar.ends.size)
                  points(location.vec[i], at[i] - (ci.offset[i] * 
                    csize)/2, pch = location.pch, col = col[(i - 
                    1L)%%length(col) + 1L], cex = location.cex)
                }
            }
            else {
                if (vertical) {
                  if (sd.x > 0 & n.vec[i] >= min.n.for.ci) 
                    errorBar(x = at[i] + ci.offset[i] * csize, 
                      y = location.vec[i], lower = ci.mat[i, 
                        "LCL"], upper = ci.mat[i, "UCL"], incr = FALSE, 
                      bar.ends = ci.bar.ends, gap = ci.bar.gap, 
                      add = TRUE, horizontal = FALSE, col = col[(i - 
                        1L)%%length(col) + 1L], lwd = ci.bar.lwd, 
                      bar.ends.size = ci.bar.ends.size)
                  points(at[i] + ci.offset[i] * csize, location.vec[i], 
                    pch = location.pch, col = col[(i - 1L)%%length(col) + 
                      1L], cex = location.cex)
                }
                else {
                  if (sd.x > 0 & n.vec[i] >= min.n.for.ci) 
                    errorBar(x = location.vec[i], y = at[i] + 
                      ci.offset[i] * csize, lower = ci.mat[i, 
                      "LCL"], upper = ci.mat[i, "UCL"], incr = FALSE, 
                      bar.ends = ci.bar.ends, gap = ci.bar.gap, 
                      add = TRUE, horizontal = TRUE, col = col[(i - 
                        1L)%%length(col) + 1L], lwd = ci.bar.lwd, 
                      bar.ends.size = ci.bar.ends.size)
                  points(location.vec[i], at[i] + ci.offset[i] * 
                    csize, pch = location.pch, col = col[(i - 
                    1L)%%length(col) + 1L], cex = location.cex)
                }
            }
        }
    }
    return.mat <- cbind(N = n.vec, Mean = location.vec, SD = scale.vec, 
        ci.mat)
    if (ci.and.test == "nonparametric") 
        dimnames(return.mat)[[2]][2:3] <- c("Median", "IQR")
    dimnames(return.mat)[[1]] <- group.names
    ret.list <- list(group.centers = at, group.stats = return.mat)
    if (p.value || (n == 2L && plot.diff)) {
        if (n == 1L) {
            p.val <- NA
            p.val.string <- ""
        }
        else if (any(sapply(groups, sd) > 0)) {
            if (n == 2) {
                if (ci.and.test == "nonparametric") {
                  test.fcn <- "wilcox.test"
                  p.val.string <- "Wilcoxon p-value"
                  if (is.null(test.arg.list)) {
                    test.arg.list <- list(conf.int = TRUE, paired = paired)
                    if (paired) 
                      p.val.string <- "Paired Wilcoxon p-value"
                  }
                  else {
                    if (all(is.na(pmatch(names(test.arg.list), 
                      "conf.int")))) {
                      test.arg.list <- c(test.arg.list, list(conf.int = TRUE))
                    }
                    else {
                      index <- (1:length(test.arg.list))[!is.na(pmatch(names(test.arg.list), 
                        "conf.int"))]
                      if (!unlist(test.arg.list[index])) 
                        test.arg.list[[index]] <- TRUE
                    }
                    if (!all(is.na(pmatch(names(test.arg.list), 
                      "paired")))) {
                      index <- (1:length(test.arg.list))[!is.na(pmatch(names(test.arg.list), 
                        "paired"))]
                      test.arg.list[[index]] <- paired
                      if (unlist(test.arg.list[index])) 
                        p.val.string <- "Paired Wilcoxon p-value"
                    }
                  }
                }
                else {
                  test.fcn <- "t.test"
                  p.val.string <- "t-test p-value"
                  if (is.null(test.arg.list)) {
                    test.arg.list <- list(var.equal = TRUE, paired = paired)
                    if (paired) 
                      p.val.string <- "Paired t-test p-value"
                  }
                  else {
                    if (!all(is.na(pmatch(names(test.arg.list), 
                      "paired")))) {
                      index <- (1:length(test.arg.list))[!is.na(pmatch(names(test.arg.list), 
                        "paired"))]
                      test.arg.list[[index]] <- paired
                      if (unlist(test.arg.list[index])) 
                        p.val.string <- "Paired t-test p-value"
                    }
                    else {
                      if (is.null(test.arg.list) || all(is.na(pmatch(names(test.arg.list), 
                        "var.equal")))) {
                        test.arg.list <- c(test.arg.list, list(var.equal = TRUE))
                      }
                      else {
                        index <- (1:length(test.arg.list))[!is.na(pmatch(names(test.arg.list), 
                          "var.equal"))]
                        if (!unlist(test.arg.list[index])) 
                          p.val.string <- "Welch t-test p-value"
                      }
                    }
                  }
                }
                alternative <- match.arg(alternative, c("two.sided", 
                  "less", "greater"))
                if (plot.diff && alternative != "two.sided") 
                  stop(paste("When there are 2 groups and plot.diff=TRUE,", 
                    "you must set alternative='two.sided'"))
                if (alternative != "two.sided") 
                  p.val.string <- paste(p.val.string, " (alternative='", 
                    alternative, "')", sep = "")
                test.list <- do.call(test.fcn, args = c(list(x = groups[[2]], 
                  y = groups[[1]], alternative = alternative, 
                  conf.level = group.difference.conf.level), 
                  test.arg.list))
                p.val <- test.list$p.value
                ci <- test.list$conf.int
                if (ci.and.test == "parametric") {
                  if (!paired) {
                    loc.diff <- -diff(test.list$estimate)
                    names(loc.diff) <- paste(paste("mean", rev(group.names)), 
                      collapse = " - ")
                  }
                  else {
                    loc.diff <- test.list$estimate
                    names(loc.diff) <- paste("mean of (", paste(rev(group.names), 
                      collapse = " - "), ") paired differences", 
                      sep = "")
                  }
                }
                else {
                  loc.diff <- test.list$estimate
                }
                if (plot.diff) {
                  usr <- par("usr")
                  if (!paired) {
                    if (is.null(diff.lim)) 
                      diff.lim <- range(pretty(range(ci)))
                    if (is.null(diff.axis.label)) 
                      diff.axis.label <- "Difference Between Groups"
                    if (vertical) {
                      if (sep.line) {
                        if (show.ci) 
                          v <- mean(c(at[2] + ci.offset[2] * 
                            csize, at[3]))
                        else v <- mean(at[2:3])
                        abline(v = v, lty = sep.lty, lwd = sep.lwd, 
                          col = sep.col)
                      }
                      par(usr = c(usr[1:2], diff.lim))
                      axis(4, at = diff.at)
                      mtext(diff.axis.label, side = 4, line = 3)
                      errorBar(x = at[3], y = loc.diff, lower = ci[1], 
                        upper = ci[2], incr = FALSE, bar.ends = ci.bar.ends, 
                        gap = ci.bar.gap, add = TRUE, horizontal = FALSE, 
                        col = diff.col, lwd = ci.bar.lwd, bar.ends.size = ci.bar.ends.size)
                      points(at[3], loc.diff, pch = location.pch, 
                        col = diff.col, cex = location.cex)
                    }
                    else {
                      if (sep.line) {
                        if (show.ci) 
                          h <- mean(c(at[2] + ci.offset[2] * 
                            csize, at[3]))
                        else h <- mean(at[2:3])
                        abline(h = h, lty = sep.lty, lwd = sep.lwd, 
                          col = sep.col)
                      }
                      par(usr = c(diff.lim, usr[3:4]))
                      axis(3, at = diff.at)
                      mtext(diff.axis.label, side = 3, line = 3)
                      errorBar(y = at[3], x = loc.diff, lower = ci[1], 
                        upper = ci[2], incr = FALSE, bar.ends = ci.bar.ends, 
                        gap = ci.bar.gap, add = TRUE, horizontal = TRUE, 
                        col = diff.col, lwd = ci.bar.lwd, bar.ends.size = ci.bar.ends.size)
                      points(loc.diff, at[3], pch = location.pch, 
                        col = diff.col, cex = location.cex)
                    }
                  }
                  else {
                    y.at <- at[3] - (ci.offset[3] * csize)/2
                    x <- groups[[2]] - groups[[1]]
                    y <- rep.int(y.at, length(x))
                    if (method == 2L) 
                      y <- y + stats::runif(length(y), -jitter, 
                        jitter)
                    else if (method == 3L) {
                      xg <- split(x, factor(x))
                      xo <- lapply(xg, seq_along)
                      x <- unlist(xg, use.names = FALSE)
                      y <- rep.int(y.at, length(x)) + (unlist(xo, 
                        use.names = FALSE) - 1) * offset * csize
                    }
                    if (is.null(diff.lim)) 
                      diff.lim <- range(pretty(c(x, range(ci))), 
                        na.rm = TRUE)
                    if (is.null(diff.axis.label)) 
                      diff.axis.label <- "Paired Difference"
                    if (vertical) {
                      if (sep.line) {
                        if (show.ci) 
                          v <- mean(c(at[2] + ci.offset[2] * 
                            csize, at[3]))
                        else v <- mean(at[2:3])
                        abline(v = v, lty = sep.lty, lwd = sep.lwd, 
                          col = sep.col)
                      }
                      par(usr = c(usr[1:2], diff.lim))
                      axis(4, at = diff.at)
                      mtext(diff.axis.label, side = 4, line = 3)
                      points(y, x, col = diff.col, pch = diff.pch, 
                        cex = points.cex, ...)
                      x <- x[is.finite(x)]
                      n.x <- sum(is.finite(x))
                      if (n.text != "none") {
                        side <- ifelse(n.text == "top", 3, 1)
                        mtext(paste("n=", n.x, sep = ""), side = side, 
                          line = n.text.line, at = at[3], cex = n.text.cex)
                      }
                      if (n.x > 0) {
                        loc <- do.call(location.fcn, args = list(x = x))
                        sca <- do.call(scale.fcn, args = list(x = x))
                        if (location.scale.text != "none") {
                          side <- ifelse(location.scale.text == 
                            "top", 3, 1)
                          dum <- format(round(c(loc, sca), location.scale.digits), 
                            nsmall = nsmall)
                          if (ci.and.test == "parametric") {
                            string1 <- "Mean="
                            string2 <- "\nSD   ="
                          }
                          else {
                            string1 <- "Median="
                            string2 <- "\nIQR     ="
                          }
                          mtext(paste(string1, dum[1], string2, 
                            dum[2], sep = ""), side = side, line = location.scale.text.line, 
                            at = at[3], cex = location.scale.text.cex)
                        }
                      }
                      if (show.ci) {
                        errorBar(x = at[3] + (ci.offset[3] * 
                          csize)/2, y = loc.diff, lower = ci[1], 
                          upper = ci[2], incr = FALSE, bar.ends = ci.bar.ends, 
                          gap = ci.bar.gap, add = TRUE, horizontal = FALSE, 
                          col = diff.col, lwd = ci.bar.lwd, bar.ends.size = ci.bar.ends.size)
                        points(at[3] + (ci.offset[3] * csize)/2, 
                          loc.diff, pch = location.pch, col = diff.col, 
                          cex = location.cex)
                      }
                    }
                    else {
                      if (sep.line) {
                        if (show.ci) 
                          h <- mean(c(at[2] + ci.offset[2] * 
                            csize, at[3]))
                        else h <- mean(at[2:3])
                        abline(h = h, lty = sep.lty, lwd = sep.lwd, 
                          col = sep.col)
                      }
                      par(usr = c(diff.lim, usr[3:4]))
                      axis(3, at = diff.at)
                      mtext(diff.axis.label, side = 3, line = 3)
                      points(x, y, col = diff.col, pch = diff.pch, 
                        cex = points.cex, ...)
                      if (show.ci) {
                        errorBar(y = at[3] + (ci.offset[3] * 
                          csize)/2, x = loc.diff, lower = ci[1], 
                          upper = ci[2], incr = FALSE, bar.ends = ci.bar.ends, 
                          gap = ci.bar.gap, add = TRUE, horizontal = TRUE, 
                          col = diff.col, lwd = ci.bar.lwd, bar.ends.size = ci.bar.ends.size)
                        points(loc.diff, at[3] + (ci.offset[3] * 
                          csize)/2, pch = location.pch, col = diff.col, 
                          cex = location.cex)
                      }
                    }
                  }
                }
            }
            else {
                y <- unlist(groups)
                group.n <- sapply(groups, length)
                group <- factor(rep(1:n, times = group.n))
                if (ci.and.test == "parametric") {
                  p.val.string <- "Anova p-value"
                  dum.aov <- aov(y ~ group)
                  dum.mat <- unclass(summary(dum.aov))[[1]]
                  p.val <- dum.mat[1, "Pr(>F)"]
                }
                else {
                  p.val.string <- "Kruskal-Wallis p-value"
                  test.list <- kruskal.test(y ~ group)
                  p.val <- test.list$p.value
                }
            }
            if (p.value & !is.na(p.val)) {
                p.val.to.show <- round(p.val, p.value.digits)
                p.val.to.show <- ifelse(p.val.to.show == 0, paste("<", 
                  format(5 * 10^-(p.value.digits + 1), scientific = FALSE)), 
                  paste("=", p.val.to.show))
                string <- paste(p.val.string, p.val.to.show)
                if (n == 2 & group.difference.ci) {
                  string1 <- ifelse(ci.and.test == "parametric", 
                    "% CI for Difference in Means: [", "% CI for Difference in Locations: [")
                  string <- paste(string, ";  ", round(100 * 
                    group.difference.conf.level, 0), string1, 
                    round(ci[1], group.difference.digits), ", ", 
                    round(ci[2], group.difference.digits), "]", 
                    sep = "")
                }
                mtext(string, line = p.value.line, cex = p.value.cex)
            }
        }
        else {
            p.val.string <- "p-value"
            p.val <- NA
            ci <- c(LCL = NA, UCL = NA)
            if (p.value) 
                warning("Constant values within each group, so between-group test not possible")
        }
        attr(p.val, "type") <- p.val.string
        if (n == 2) {
            ret.list <- c(ret.list, list(group.difference = loc.diff))
        }
        ret.list <- c(ret.list, list(group.difference.p.value = p.val))
        if (n == 2) {
            ret.list <- c(ret.list, list(group.difference.conf.int = ci))
        }
    }
    invisible(ret.list)
}
