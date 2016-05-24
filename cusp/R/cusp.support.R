`cusp.support` <-
function (alpha, beta, eps = 1e-05) 
unlist(sapply(1:2, function(i) optimize(function(x) {
    v = dcusp.unnorm(x, alpha, beta)
    v + (v < abs(eps)) * x^2
}, low = c(-999, cusp.extrema(alpha, beta)[3])[i], upp = c(cusp.extrema(alpha, 
    beta)[1], 999)[i]))[1, ])

