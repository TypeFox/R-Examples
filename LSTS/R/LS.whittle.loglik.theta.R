LS.whittle.loglik.theta = function(x, series, order = c(p = 0, q = 0), ar.order = NULL, ma.order = NULL, sd.order = NULL, d.order = NULL, include.d = FALSE, N = NULL, S = NULL, include.taper = TRUE,  sd.par = 1){
  x = c(x,sd.par)
  LS.whittle.loglik(x = x, series = series, order = order, ar.order = ar.order, ma.order = ma.order, sd.order = sd.order, d.order = d.order, include.d = include.d, N = N, S = S, include.taper = include.taper)
}
