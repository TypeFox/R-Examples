`ARMAspec` <-
function (model, freq = seq(0, 0.5, 0.001), plot = TRUE, ...) 
{
    comp.spec = function(p, ar, freq, period = 1) {
        arg1 = outer(freq, 1:p, function(x, y) {
            2 * pi * x * y * period
        })
        cosy = cbind(1, cos(arg1)) %*% c(1, ar)
        siny = sin(arg1) %*% ar
        cosy * cosy + siny * siny
    }
    spec = freq * 0 + 1
    if (!is.list(model)) 
        stop("'model' must be list")
    p <- length(model$ar)
    if (p) {
        minroots <- min(Mod(polyroot(c(1, -model$ar))))
        if (minroots <= 1) 
            stop("'ar' part of model is not stationary")
        spec = spec/comp.spec(p, -model$ar, freq)
    }
    q <- length(model$ma)
    if (q) {
        spec = spec * comp.spec(q, model$ma, freq)
    }
    if (!is.null(model$seasonal)) {
        seasonal = model$seasonal
        if (!is.null(seasonal$period)) 
            period = seasonal$period
        else period = 12
        P = length(seasonal$sar)
        if (P) {
            minroots <- min(Mod(polyroot(c(1, -seasonal$sar))))
            if (minroots <= 1) 
                stop("'sar' part of model is not stationary")
            spec = spec/comp.spec(P, -seasonal$sar, freq, period)
        }
        Q = length(seasonal$sma)
        if (Q) {
            spec = spec * comp.spec(Q, seasonal$sma, freq, period)
        }
    }
    if (!is.null(model$sigma2)) 
        spec = spec * model$sigma2
    ylim = c(0, max(spec))
    if (plot) {
        plot(y = spec, x = freq, ylim = ylim, xlab = "Frequency", 
            ylab = "Spectral Density", type = "l", ...)
    abline(h = 0)}
    invisible(list(spec = spec, freq = freq, model = model))
}

