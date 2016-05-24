`aweibull7` <-
function(p, lower = 0, upper = 365) {
     integrate(f = fweibull7, lower = lower, upper = upper, p = p)$value
}

