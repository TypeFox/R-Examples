`plot.Arima` <-
function (x, n.ahead = 12, col = "black", ylab = object$series, 
    lty = 2, n1, newxreg, transform, Plot = TRUE, ...) 
{
    object = x
    if (missing(newxreg)) 
        pred.res = predict(object, n.ahead = n.ahead)
    else pred.res = predict(object, n.ahead = n.ahead, newxreg = newxreg)
    pred = pred.res$pred
    lpi = pred - 1.96 * pred.res$se
    upi = pred + 1.96 * pred.res$se
    x = eval(object$call$x)
    if (missing(n1)) 
        n1 = start(x)
    n2 = end(x)
    x = ts(window(x, start = n1, end = n2), start = n1, end = n2, 
        frequency = frequency(x))
    y = ts(c(x, pred), start = start(x), frequency = frequency(x))
    z = ts(window(y, start = n2), start = n2, frequency = frequency(x))
    if (!missing(transform)) {
        y = transform(y)
        x = transform(x)
        z = transform(z)
        upi = transform(upi)
        lpi = transform(lpi)
    }
    if (Plot) {
        mc = match.call(expand.dots = FALSE)
        ellipsis = mc$...
        if (!is.null(ellipsis$type)) 
            oldtype = ellipsis$type
        else oldtype = "b"
        ellipsis$type = "n"
        do.call(plot, c(list(x = y, ylim = range(c(y, lpi, upi)), 
            ylab = ylab), ellipsis))
        ellipsis$type = oldtype
        if (!is.null(ellipsis$pch)) 
            old.pch = ellipsis$pch
        else old.pch = 1
        ellipsis$pch = NA
        do.call(lines, c(list(x = z, lty = lty), ellipsis))
        ellipsis$pch = 1
        ellipsis$type = "o"
        do.call(lines, c(list(x = x, lty = 1), ellipsis))
        ellipsis$type = oldtype
        ellipsis$pch = old.pch
        do.call(points, c(list(x = window(y, start = n2 + c(0, 
            1 * deltat(y))), lty = lty), ellipsis))
        lines(upi, lty = 3, col = col)
        lines(lpi, lty = 3, col = col)
    }
    invisible(list(pred = window(y, start = n2 + c(0, 1 * deltat(y))), 
        lpi = lpi, upi = upi))
}

