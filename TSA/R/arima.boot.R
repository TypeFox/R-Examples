arima.boot <-
function (arima.fit, cond.boot = FALSE, is.normal = TRUE, B = 1000, 
    init, ntrans = 100) 
{
    arimab.sim <- function(model = list(ar = NA, ma = NA, d = 0), 
        n, ntrans = 100, init = rep(0, 100), cond.boot = FALSE, 
        residuals, include.mean = TRUE, mean, is.normal = TRUE, 
        sd, ...) {
        if (!is.null(model$ar)) 
            ar = model$ar
        else ar = NA
        if (!is.null(model$ma)) 
            ma = model$ma
        else ma = NA
        if (!is.null(model$d)) 
            d = model$d
        else d = 0
        if (include.mean) 
            intercept = mean
        else intercept = 0
        if (!any(is.na(ar))) {
            d1 = d
            if (d1 > 0) {
                bar = c(1, -ar)
                while (d1 > 0) {
                  bar = c(bar, 0) - c(0, bar)
                  d1 = d1 - 1
                }
                ar = -bar[-1]
            }
            p = length(ar)
        }
        else {
            p = 0
            d1 = d
            if (d1 > 0) {
                ar = c(1, 0)
                while (d1 > 0) {
                  bar = c(bar, 0) - c(0, bar)
                  d1 = d1 - 1
                }
                ar = -bar[-1]
            }
            if (!any(is.na(ar))) 
                p = length(ar)
        }
        if (!is.na(ma)) 
            q = length(ma)
        else q = 0
        initial = rev(init[1:p]) - intercept
        if (cond.boot) {
            ntrans = 0
        }
        if (is.normal) {
            if (q > 0) 
                noise = filter(rnorm(n = n + ntrans, mean = 0, 
                  sd = sd), init = rnorm(n = q, mean = 0, sd = sd), 
                  filter = ma, method = "convolution", sides = 1)
            else noise = rnorm(n = n + ntrans, mean = 0, sd = sd)
        }
        else {
            if (q > 0) 
                noise = filter(sample(residuals, replace = TRUE, 
                  size = n + ntrans), init = sample(residuals, 
                  size = q, replace = TRUE), filter = ma, method = "convolution", 
                  sides = 1)
            else noise = sample(residuals, size = n + ntrans, 
                replace = TRUE)
        }
        boot = filter(noise, filter = ar, method = "recursive", 
            init = initial, sides = 1) + intercept
        boot = boot[(ntrans + 1):(n + ntrans)]
        if (cond.boot) 
            boot = c(rev(initial) + intercept, rev(rev(boot)[-seq(initial)]))
        boot
    }

    if (!missing(arima.fit)) {
        order = arima.fit$call$order
        p = order[[2]]
        d = order[[3]]
        q = order[[4]]
        order = c(p, d, q)
        n = length(eval(arima.fit$call[[2]]))
        if (p > 0) 
            ar = coef(arima.fit)[1:p]
        else ar = NA
        if (q > 0) 
            ma = coef(arima.fit)[(p + 1):(p + q)]
        else ma = NA
        model = list(ar = ar, ma = ma, d = d)
        if (any(is.null(arima.fit$call$include.mean))) 
            include.mean = TRUE
        else include.mean = (arima.fit$call$include.mean) == 
            "T"
        if (include.mean) 
            mean = coef(arima.fit)[(p + q + 1)]
        else mean = 0
        sd = arima.fit$sigma2^0.5
        residuals = residuals(arima.fit)
        if (missing(init)) {
            init = eval(arima.fit$call[[2]])[1:(p + d)]
        }
    }
    else stop("Need to input the fitted model")
    coefm = NULL
    i=1
    repeat {
        boot = arimab.sim(model = model, n = n, ntrans = ntrans, 
            init = init, cond.boot = cond.boot, residuals = residuals, 
            include.mean = include.mean, mean = mean, is.normal = is.normal, 
            sd = sd)
        res = try(arima(boot, order = order, include.mean = include.mean), silent=TRUE)
        if(class(res)[1] == "try-error") next        
        if (any(is.na(res))) next 
            coefm = rbind(coefm, c(coef(res), res$sigma2))
            i=i+1
            if(i>B) break
    }
    colnames(coefm) = c(names(coef(res)), "noise var")
    invisible(coefm)
}
