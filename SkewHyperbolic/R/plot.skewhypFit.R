######### Plot method for skewhypFit ########################################
plot.skewhypFit <- function(x, which = 1:4,
                            plotTitles = paste(c("Histogram of ",
                            "Log-Histogram of ", "Q-Q Plot of ",
                            "P-P Plot of "), x$obsName, sep = ""),
                            ask = prod(par("mfcol")) < length(which) &&
                            dev.interactive(), ...) {

    if (!class(x) == "skewhypFit")
        stop("Object must belong to class skewhypFit")

    if (ask) {
        op <- par(ask = TRUE)
        on.exit(par(op))
    }

    par(mar = c(6, 4, 4, 2) + 0.1)
    show <- rep(FALSE, 4)
    show[which] <- TRUE
    param <- x$param
    breaks <- x$breaks
    empDens <- x$empDens
    mipoints <- x$midpoints
    obs <- x$x
    obsName <- x$xName
    skewhypDens <- function(obs) {
        dskewhyp(obs, param = param, tolerance = .Machine$double.eps ^ 0.5)
    }
    logskewhypDens <- function(obs) {
        dskewhyp(obs, param = param, tolerance = .Machine$double.eps ^ 0.5,
                 log = TRUE)
    }
    ymax <- 1.06 * max(skewhypDens(seq(min(breaks), max(breaks),
                                       length = length(obs))),
                       empDens, na.rm = TRUE)

    if (show[1]) {#histogram
        hist(obs, breaks = breaks, freq = FALSE, ylim = c(0,ymax),
             main = paste(plotTitles[1], obsName), ...)
        curve(skewhypDens(x), add = TRUE, col = 2,  ...)
        title(sub = paste("param = (", round(param[1], 3), ",",
              round(param[2], 3), ",", round(param[3], 3), ",",
              round(param[4], 3), ")", sep = ""), ...)
    }
    if (show[2]) {#log histogram
        logHist(obs, breaks, include.lowest = TRUE, right = FALSE,
                main = paste(plotTitles[2], obsName), ...)
        curve(logskewhypDens, min(breaks) - 1, max(breaks) + 1, add = TRUE,
              ylab = NULL, xlab = NULL, col = 2, ...)
        title(sub = paste("param = (", round(param[1], 3), ",",
              round(param[2], 3), ",", round(param[3], 3), ",",
              round(param[4], 3), ")", sep = ""), ...)
    }
    if (show[3]) {#qq
        qqskewhyp(obs, param = param, main = paste(plotTitles[3], obsName),
                  line = TRUE, ...)
    }
    if (show[4]) {#pp
        ppskewhyp(obs, param = param, main = paste(plotTitles[4], obsName),
                  line = TRUE, ...)
    }
    invisible()
}
###### Print method for skewhypFit #####################################
print.skewhypFit <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    if (!class(x) == "skewhypFit")
        stop("Object must belong to class skewhypFit")

    cat("\nData:     ", x$xName, "\n")
    cat("Parameter estimates:\n")

    if (is.null(x$sds)) {
        print.default(format(x$param, digits = digits),
            print.gap = 2, quote = FALSE)
    }else {
        ans <- format(rbind(x$param, x$sds), digits = digits)
        ans[1, ] <- sapply(ans[1, ], function(obs) paste("", obs))
        ans[2, ] <- sapply(ans[2, ], function(obs) paste("(", obs,
            ")", sep = ""))
        dn <- dimnames(ans)
        dn[[1]] <- rep("", 2)
        dn[[2]] <- paste(substring("      ", 1, (nchar(ans[2,
            ]) - nchar(dn[[2]]))%/%2), dn[[2]])
        dn[[2]] <- paste(dn[[2]], substring("      ", 1, (nchar(ans[2,
            ]) - nchar(dn[[2]]))%/%2))
        dimnames(ans) <- dn
        print.default(ans, print.gap = 2, quote = FALSE)
    }

    cat("Likelihood:        ", x$maxLik, "\n")
    cat("Method:            ", x$method, "\n")
    cat("Convergence code:  ", x$conv, "\n")
    cat("Iterations:        ", x$iter, "\n")
    invisible(x)
}
