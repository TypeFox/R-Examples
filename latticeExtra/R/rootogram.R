



prepanel.rootogram <-
    function(x, y = table(x),
             dfun = NULL,
             transformation = sqrt,
             hang = TRUE,
             probability = TRUE,
             ...)
{
    stopifnot(is.function(dfun))
    if (probability) y <- y / sum(y)
    yy <- transformation(y)
    xx <- sort(unique(x))
    dotArgs <- list(...)
    dfunArgs <- names(formals(dfun))
    if (!("..." %in% dfunArgs))
        dotArgs <- dotArgs[dfunArgs[-1]]
    dd <- transformation(do.call(dfun, c(list(xx), dotArgs)))
    list(xlim = range(xx),
         ylim =
         if (hang) range(dd, dd-yy, 0)
         else range(dd, yy, 0),
         dx = diff(xx),
         dy = diff(dd))
}


panel.rootogram <-
    function(x, y = table(x),
             dfun = NULL,
             col = plot.line$col,
             lty = plot.line$lty,
             lwd = plot.line$lwd,
             alpha = plot.line$alpha,
             transformation = sqrt,
             hang = TRUE,
             probability = TRUE,
             type = "l", pch = 16,
             ...)
{
    plot.line <- trellis.par.get("plot.line")
    ref.line <- trellis.par.get("reference.line")
    stopifnot(is.function(dfun))
    if (probability) y <- y / sum(y)
    yy <- transformation(y)
    xx <- sort(unique(x))
    dotArgs <- list(...)
    dfunArgs <- names(formals(dfun))
    if (!("..." %in% dfunArgs))
        dotArgs <- dotArgs[dfunArgs[-1]]
    dd <- transformation(do.call(dfun, c(list(xx), dotArgs)))
    panel.abline(h = 0,
                 col = ref.line$col,
                 lty = ref.line$lty,
                 lwd = ref.line$lwd,
                 alpha = ref.line$alpha)
    panel.segments(xx,
                   if (hang) dd else 0,
                   xx,
                   if (hang) (dd - yy) else yy,
                   col = col,
                   lty = lty,
                   lwd = lwd,
                   alpha = alpha,
                   ...)
    if ("l" %in% type) panel.lines(xx, dd)
    if ("p" %in% type) panel.points(xx, dd, pch = pch)
}


rootogram <-
    function(x, ...)
    UseMethod("rootogram")





rootogram.formula <-
    function(x, data = parent.frame(),
             ylab = expression(sqrt(P(X == x))),
             prepanel = prepanel.rootogram,
             panel = panel.rootogram,
             ...,
             probability = TRUE)
{
    if (!probability && missing(ylab)) ylab <- NULL
    if (length(x) == 2) ## formula like ~ x
        foo <-
            densityplot(x, data,
                        prepanel = prepanel,
                        panel = panel,
                        ylab = ylab,
                        ...,
                        probability = probability)
    else ## formula like y ~ x 
        foo <-
            xyplot(x, data,
                   prepanel = prepanel,
                   panel = panel,
                   ylab = ylab,
                   ...,
                   probability = probability)
    foo$call <- sys.call(sys.parent()); foo$call[[1]] <- quote(rootogram)
    foo
}

