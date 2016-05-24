my.approx <-
function (x, y = NULL, xout, rule = 3, ...) 
{
    if (!is.null(y) && rule == 3 && (any(xout < min(x)) | any(xout > 
        max(x)))) {
        if (length(x) != length(y)) 
            warning("x and y length differ")
        xl <- xout[xout < min(x)]
        xm <- xout[min(x) <= xout & xout <= max(x)]
        xr <- xout[xout > max(x)]
        beta1 = (max(y) - min(y))/(max(x) - min(x))
        beta0 = max(y) - beta1 * max(x)
        yl <- vector(length = length(xl))
        yr <- vector(length = length(xr))
        for (i in seq(along = xl)) yl[i] <- beta1 * xl[i] + beta0
        for (i in seq(along = xr)) yr[i] <- beta1 * xr[i] + beta0
        Erg <- approx(x = x, y = y, xout = xm, method = "linear")
        list(x = c(xl, Erg$x, xr), y = c(yl, Erg$y, yr))
    }
    else {
        approx(x, y, xout, rule = rule, ...)
    }
}
