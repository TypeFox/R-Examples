plot.sbchoice <- function (x, type = NULL, main = NULL, sub = NULL,
    xlab = NULL, ylab = NULL, lwd = NULL, lty = NULL, xlim = NULL,
    ylim = NULL, bid = NULL, ...)
{

    if (is.null(bid)) {
        minbid <- min(x$bid)
        maxbid <- max(x$bid)
    } else {
        minbid <- min(bid)
        maxbid <- max(bid)
    }
    bidseq <- seq(minbid, maxbid, by = (maxbid - minbid)/100)
    b <- x$coefficients
    npar <- length(b)
    Xcov <- colMeans(x$covariates)
    Vcov <- sum(Xcov * b[-npar])
    V <- Vcov + bidseq * b[npar]
    dist <- x$distribution

    if(dist == "logistic" | dist == "log-logistic") {
        pr <- plogis(-V, lower.tail = FALSE, log.p = FALSE)
    } else if (dist == "normal" | dist == "log-normal") {
        pr <- pnorm(-V, lower.tail = FALSE, log.p = FALSE)
    } else if (dist == "weibull") {
        pr <- pweibull(exp(-V), shape = 1, lower.tail = FALSE, log.p = FALSE)
    }

    if (is.null(type)) type <- "l"
    if (is.null(main)) main <- ""
    if (is.null(sub))  sub <- ""
    if (is.null(xlab)) xlab <- names(b[npar])
    if (is.null(ylab)) ylab <- "Probability of selecting yes"
    if (is.null(lwd))  lwd <- 3
    if (is.null(lty))  lty <- 1
    if (is.null(xlim)) xlim <- c(0.96 * minbid, 1.04 * maxbid)
    if (is.null(ylim)) ylim <- c(0, 1)

    plot(x = bidseq, y = pr, xlab = xlab, ylab = ylab, main = main,
         lwd = lwd, type = type, xlim = xlim, ylim = ylim, ...)

    invisible(list(x = bidseq, utility = V, probability = pr))
}
