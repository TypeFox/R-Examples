tsdiag.Arima <-
function (object, gof.lag, tol = 0.1, col = "red", omit.initial = TRUE, 
    ...) 
{
    opar=par(mfrow = c(3, 1), mar = c(3, 4, 3, 2) + 0.1, oma = c(1, 
        0, 2, 0))
    n = length(eval(object$call[[2]]))
    if (missing(gof.lag)) 
        lag.max = 10 * log10(n)
    else lag.max = gof.lag
    phi = object$mod$phi
    Delta = c(1, -object$mod$Delta)
    theta = object$mod$theta
    p1 = length(phi)
    d1 = length(Delta)
    q = length(theta)
    residuals = residuals(object)
    if (omit.initial) 
        residuals = window(residuals, start = time(residuals)[d1 + 
            1])
    std.res = residuals/object$sigma2^0.5
    n = length(std.res)
    h1 = qnorm(0.025/n)
    plot(std.res, ylab = "Standardized Residuals", type = "p", 
        ...)
    abline(h = h1, lty = 2, col = col)
    abline(h = -h1, lty = 2, col = col)
    abline(h = 0)
    acf(as.numeric(residuals), lag.max = lag.max, ylab = "ACF of Residuals", 
        ci.col = col, main = "", ...)
    ar = phi
    ma = c(1, theta)
    psiv = rep(0, lag.max)
    psiv[seq(ma)] = ma
    if (p1 > 0) 
        psiv = filter(psiv, filter = ar, method = "recursive", 
            sides = 1)
    psiv = psiv[-1]
    test = abs(psiv) < tol
    test = rev(cumprod(rev(test)))
    if (any(test==1)) 
        k = (seq(along.with=test)[test == 1])[1]
    else stop("increase gof.lag as psi weights are not small enough for the Ljung-Box tests")
    lbv = rep(NA, lag.max)
    for (i in k:lag.max) {
        lbv[i] = LB.test(object, lag = i, no.error = TRUE, omit.initial = omit.initial)$p.value
    }
    plot(y = lbv, x = 1:lag.max, ylim = c(0, 1), pch = 21, ylab = "P-values", 
        xlab = "Number of lags", ...)
    abline(h = 0.05, lty = 2, col = col)
    par(opar)
    invisible()
}