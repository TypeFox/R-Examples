`cusp.nc` <-
Vectorize(function (alpha, beta, mom.order = 0, ...) 
tryCatch(integrate(function(x) x^mom.order * exp(alpha * x + 
    beta * x^2/2 - x^4/4), -Inf, Inf, ...)$value, error = function(e) Inf))

