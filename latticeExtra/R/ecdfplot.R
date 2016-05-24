
prepanel.ecdfplot <-
    function(x, f.value = NULL, ...)
{
    ans <-
        prepanel.default.qqmath(x,
                                f.value = f.value,
                                distribution = qunif)
    with(ans, list(xlim = ylim, ylim = c(0, 1),
                   dx = dy, dy = dx))
}



panel.ecdfplot <-
    function(x, f.value = NULL, type = "s",
             groups = NULL, qtype = 7,
             ref = TRUE,
             ...)
{
    if (ref)
    {
        reference.line <- trellis.par.get("reference.line")
        do.call(panel.abline, c(list(h = c(0, 1)), reference.line))
    }
    x <- as.numeric(x)
    distribution <- qunif
    nobs <- sum(!is.na(x))
    if (!is.null(groups))
    {
        panel.superpose(x, y = NULL,
                        f.value = f.value, type = type,
                        distribution = distribution, 
                        qtype = qtype,
                        groups = groups,
                        panel.groups = panel.ecdfplot,
                        ...)
    }
    else if (nobs)
    {
        if (is.null(f.value))
        {
            panel.xyplot(x = sort(x),
                         y = seq_len(nobs) / nobs,
                         type = type, 
                         ...)
        }
        else
        {
            p <- if (is.numeric(f.value)) f.value else f.value(nobs)
            panel.xyplot(x = quantile(x, p, names = FALSE, type = qtype, na.rm = TRUE),
                         y = distribution(p),
                         type = type, 
                         ...)
        }
    }
}



ecdfplot <- 
    function (x, data, ...) 
    UseMethod("ecdfplot")


ecdfplot.formula <- 
    function (x, data = NULL,
              prepanel = "prepanel.ecdfplot", 
              panel = "panel.ecdfplot",
              ylab = gettext("Empirical CDF"),
              ...) 
{
    ccall <- match.call()
    ocall <- sys.call(sys.parent()); ocall[[1]] <- quote(ecdfplot) ## for nice $call
    ccall$data <- data
    ccall$prepanel <- prepanel
    ccall$panel <- panel
    ccall$ylab <- ylab
    ccall[[1]] <- quote(lattice::densityplot)
    ans <- eval.parent(ccall)
    ans$call <- ocall
    ans
}


ecdfplot.numeric <- 
    function (x, data = NULL, xlab = deparse(substitute(x)), ...)
{
    ccall <- match.call()
    ocall <- sys.call(sys.parent()); ocall[[1]] <- quote(ecdfplot) ## for nice $call
    if (!is.null(ccall$data))
        warning("explicit 'data' specification ignored")
    ccall$data <- list(x = x)
    ccall$xlab <- xlab
    ccall$x <- ~x
    ccall[[1]] <- quote(latticeExtra::ecdfplot)
    ans <- eval.parent(ccall)
    ans$call <- ocall
    ans
}

