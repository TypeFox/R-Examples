### Adopted from plot.spec R code
 spec.ci <- function (spec.obj, coverage = 0.95)
    {
        ## A utility function for plot.spec which calculates the confidence
        ## interval (centred around zero). We use a conditional argument to
        ## ensure that the ci always contains zero.

        if (coverage < 0 || coverage >= 1)
            stop("coverage probability out of range [0,1)")
        tail <- (1 - coverage)
        df <- spec.obj$df
        upper.quantile <- 1 - tail * pchisq(df, df, lower.tail = FALSE)
        lower.quantile <- tail * pchisq(df, df)
        1/(qchisq(c(upper.quantile, lower.quantile), df)/df)
    }

"plotSpecLs" <-
function (x, add = FALSE, ci = 0.95, log = c("yes", "dB", "no"), 
    xlab = "frequency", ylab = NULL, type = "l",
    main = NULL, sub = NULL, ...) 
{
    log <- match.arg(log)
    if (is.null(ylab)) 
        ylab <- if (log == "dB") 
            "spectrum (dB)"
        else "spectrum"
    if (is.logical(log)) 
        log <- if (log) 
            "yes"
        else "no"
    if (missing(log) && getOption("ts.S.compat")) 
        log <- "dB"
    log <- match.arg(log)
    ylog <- ""
    if (log == "dB") 
        x$spec <- 10 * log10(x$spec)
    if (log == "yes") 
        ylog <- "y"
    if (add) {
        matplot(x$freq, x$spec, type = type, add = TRUE, ...)
    }
    else {
        matplot(x$freq, x$spec, xlab = xlab, ylab = ylab, type = type, 
            log = ylog, ...)
        if (ci <= 0 || !is.numeric(x$df) || log == "no") {
            ci.text <- ""
        }
        else {
            conf.lim <- spec.ci(x, coverage = ci)
            if (log == "dB") {
                conf.lim <- 10 * log10(conf.lim)
                conf.y <- max(x$spec) - conf.lim[2]
                conf.x <- max(x$freq) - x$bandwidth
                lines(rep(conf.x, 2), conf.y + conf.lim, ...)
                lines(conf.x + c(-0.5, 0.5) * x$bandwidth, rep(conf.y, 
                  2))
                ci.text <- paste(", ", round(100 * ci, 2), "% C.I. is (", 
                  paste(format(conf.lim, digits = 3), collapse = ","), 
                  ")dB", sep = "")
            }
            else {
                ci.text <- ""
                conf.y <- max(x$spec)/conf.lim[2]
                conf.x <- max(x$freq) - x$bandwidth
                lines(rep(conf.x, 2), conf.y * conf.lim)
                lines(conf.x + c(-0.5, 0.5) * x$bandwidth, rep(conf.y, 
                  2))
            }
        }
        if (is.null(main)) 
            main <- paste(if (!is.null(x$series)) 
                paste("Series:", x$series)
            else "Lomb-Scargle peridogram", x$method, sep = "\n")
        if (is.null(sub) && is.numeric(x$bandwidth)) 
            sub <- paste("bandwidth = ", format(x$bandwidth, 
                digits = 3), ci.text, sep = "")
        title(main = main, sub = sub)
    }
    invisible(x)
}
