LS.whittle = function(series, start, order = c(p = 0, q = 0), ar.order = NULL, ma.order = NULL, sd.order = NULL, d.order = NULL, include.d = FALSE, N = NULL, S = NULL, include.taper = TRUE, control = list(), lower = -Inf, upper = Inf, m = NULL, n.ahead = 0){
  
  series = c(series, rep(NA,n.ahead))
  
  aux = nlminb(start = start, objective = LS.whittle.loglik, series = series, order = order, ar.order = ar.order, ma.order = ma.order, sd.order = sd.order, d.order = d.order, include.d = include.d, N = N, S = S, include.taper = include.taper, lower = lower, upper = upper, control = control)
  
  T. = length(series)
  loglik = -aux$objective
  npar = length(aux$par)
  aic = -2 * loglik + 2 * npar/T.
  
  aux. = LS.kalman(series = series-mean(series, na.rm = TRUE), start = aux$par, order = order, ar.order = ar.order, ma.order = ma.order, sd.order = sd.order, d.order = d.order, include.d = include.d, m = m)
  
  x = aux$par
  k = npar
  if(is.null(sd.order)){sd.order = 0}
  if(is.null(d.order)){d.order = 0}
  G1 = hessian(f = LS.whittle.loglik.theta, x0 = x[1:(k-d.order-1)], series = series, order = order, ar.order = ar.order, ma.order = ma.order, sd.order = sd.order, d.order = d.order, include.d = include.d, N = N, S = S, include.taper = include.taper, sd.par = x[(k-d.order):k])
  G2 = hessian(f = LS.whittle.loglik.sd, x0 = x[(k-d.order):k], series = series, order = order, ar.order = ar.order, ma.order = ma.order, sd.order = sd.order, d.order = d.order, include.d = include.d, N = N, S = S, include.taper = include.taper, theta.par = x[1:(k-d.order-1)])
  G = matrix(0, ncol = k, nrow = k)
  G[1:(k-d.order-1),1:(k-d.order-1)] = G1
  G[(k-d.order):k, (k-d.order):k] = G2
  G = solve(G)/T.
  
  fitted.values = aux.$fitted.values + mean(series, na.rm = TRUE)
  fitted.values[is.na(series) == 1] = NA
  pred = aux.$fitted.values + mean(series, na.rm = TRUE)
  pred[is.na(series) == 0] = NA
  se = sqrt(aux.$delta)
  se[is.na(series) == 0] = NA
  
  list(coef = aux$par, var.coef = G, loglik = loglik, aic = aic, series = series, residuals = aux.$residuals, fitted.values = fitted.values, pred = pred, se = se, model = list(order = order, ar.order = ar.order, ma.order = ma.order, sd.order = sd.order, d.order = d.order, include.d = include.d, include.taper = include.taper))
}