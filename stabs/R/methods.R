print.stabsel <- function(x, decreasing = FALSE, print.all = TRUE, ...) {

    cat("\tStability Selection")
    if (x$assumption == "none")
        cat(" without further assumptions\n")
    if (x$assumption == "unimodal")
        cat(" with unimodality assumption\n")
    if (x$assumption == "r-concave")
        cat(" with r-concavity assumption\n")
    if (length(x$selected) > 0) {
        cat("\nSelected variables:\n")
        print(x$selected)
    } else {
        cat("\nNo variables selected\n")
    }
    cat("\nSelection probabilities:\n")
    if (print.all) {
        print(sort(x$max, decreasing = decreasing))
    } else {
        print(sort(x$max[x$max > 0], decreasing = decreasing))
    }
    cat("\n---\n")
    print.stabsel_parameters(x, heading = FALSE)
    cat("\n")
    invisible(x)
}

## strip stabsel results
parameters <- function(object) {
    if (!inherits(object, "stabsel") || !inherits(object, "stabsel_parameters"))

    res <- object[c("cutoff", "q", "PFER", "specifiedPFER", "p", "B",
                    "sampling.type", "assumption")]
    class(res) <- "stabsel_parameters"
    res
}

stabsel_parameters.stabsel <- function(p, ...) {
    parameters(p)
}

print.stabsel_parameters <- function(x, heading = TRUE, ...) {
    if (heading) {
        cat("Stability selection")
        if (x$assumption == "none")
            cat(" without further assumptions\n\n")
        if (x$assumption == "unimodal")
            cat(" with unimodality assumption\n\n")
        if (x$assumption == "r-concave")
            cat(" with r-concavity assumption\n\n")
    }
    cat("Cutoff: ", x$cutoff, "; ", sep = "")
    cat("q: ", x$q, "; ", sep = "")
    if (x$sampling.type == "MB")
        cat("PFER: ", x$PFER, "\n")
    else
        cat("PFER (*): ", x$PFER,
            "\n   (*) or expected number of low selection probability variables\n")
    if (!is.null(x$specifiedPFER)) {
        cat("PFER (specified upper bound): ", x$specifiedPFER, "\n")
    } else {
        if (!!is.null(x$call) && !is.null(x$call[["PFER"]])) {
            cat("PFER (specified upper bound): ", x$call[["PFER"]], "\n")
        }
    }
    p <- NULL
    if (!is.null(x$p)) {
        p <- x$p
    } else {
        if (!is.null(x$max))
            p <- length(x$max)
    }
    if (!is.null(p))
        cat("PFER corresponds to signif. level ",
            signif(x$PFER / p, 3), " (without multiplicity adjustment)\n",
            sep = "")
    invisible(x)
}

plot.stabsel <- function(x, main = deparse(x$call), type = c("maxsel", "paths"),
                         xlab = NULL, ylab = NULL,
                         col = NULL, ymargin = 10, np = sum(x$max > 0),
                         labels = NULL, ...) {

    type <- match.arg(type)

    if (is.null(col))
        col <- hcl(h = 40, l = 50, c = x$max / max(x$max) * 490)
    if (type == "paths" && is.null(x$phat)) {
        warning("Stability paths ", sQuote("x$phat"), " are missing, ",
                "plot maximum selection frequency instead")
        type <- "maxsel"
    }
    if (type == "paths") {
        ## if par(mar) not set by user ahead of plotting
        if (all(par()[["mar"]] == c(5, 4, 4, 2) + 0.1))
            ..old.par <- par(mar = c(5, 4, 4, ymargin) + 0.1)
        h <- x$phat
        h <- h[rowSums(h) > 0, , drop = FALSE]

        if (is.null(xlab)) {
            if (inherits(x, "stabsel_mboost")) {
                xlab <- "Number of boosting iterations"
            } else {
                xlab <- "Step"
            }
        }
        if (is.null(ylab)) {
            ylab <- "Selection probability"
        }
        matplot(t(h), type = "l", lty = 1,
                xlab = xlab, ylab = ylab,
                main = main, col = col[x$max > 0], ylim = c(0, 1), ...)
        abline(h = x$cutoff, lty = 1, col = "lightgray")
        if (is.null(labels))
            labels <- rownames(x$phat)
        axis(4, at = x$phat[rowSums(x$phat) > 0, ncol(x$phat)],
             labels = labels[rowSums(x$phat) > 0], las = 1)
    } else {
        ## if par(mar) not set by user ahead of plotting
        if (all(par()[["mar"]] == c(5, 4, 4, 2) + 0.1))
            ..old.par <- par(mar = c(5, ymargin, 4, 2) + 0.1)
        if (np > length(x$max))
            stop(sQuote("np"), "is set too large")
        inc_freq <- x$max  ## inclusion frequency
        if (is.null(xlab)) {
            xlab <- expression(hat(pi))
        }
        if (is.null(ylab)) {
            ylab <- ""
        }
        plot(tail(sort(inc_freq), np), 1:np,
             type = "n", yaxt = "n", xlim = c(0, 1),
             ylab = ylab, xlab = xlab,
             main = main, ...)
        abline(h = 1:np, lty = "dotted", col = "grey")
        points(tail(sort(inc_freq), np), 1:np, pch = 19,
               col = col[tail(order(inc_freq), np)])
        if (is.null(labels))
            labels <- names(x$max)
        axis(2, at = 1:np, labels[tail(order(inc_freq), np)], las = 2)
        ## add cutoff
        abline(v = x$cutoff, col = "grey")
    }
    if (exists("..old.par"))
        par(..old.par) # reset plotting settings
}

selected <- function(object, ...)
    UseMethod("selected", object)

selected.stabsel <- function(object, ...)
    object$selected
