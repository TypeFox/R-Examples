extract.ew.density <-
function(initial.values = c(NA, NA, NA), r, y, te, s0, market.calls, call.strikes, call.weights = 1, lambda = 1, hessian.flag = F, cl = list(maxit=10000))
{

  if ( sum(is.na(initial.values)) >= 1 )  
  {
    sigma.range = seq(from = 0.1, to = 0.9, by = 0.1)
    v           =  sqrt(exp(sigma.range^2 * te) - 1)
    skew.ln     =  3 * v + v^3
    kurt.ln     =  16 * v^2 + 15 * v^4 + 6 * v^6 + v^8
    ew.grid     =  cbind(sigma.range, skew.ln, kurt.ln)

    ew.vals     = apply(ew.grid, 1, ew.objective, r = r, y = y, te = te, s0 = s0, market.calls = market.calls, call.strikes = call.strikes, call.weights = call.weights, lambda = 1)
    initial.values = as.numeric(ew.grid[which.min(ew.vals),])
  }

  optim.obj = optim(initial.values, ew.objective, r=r, y = y, te=te, s0 = s0, market.calls = market.calls, 
               call.strikes = call.strikes, lambda = lambda, hessian = hessian.flag , call.weights = call.weights, control = cl )

  sigma    = optim.obj$par[1]
  skew     = optim.obj$par[2]
  kurt     = optim.obj$par[3]

  converge.result = optim.obj$convergence
  if (hessian.flag) h = optim.obj$hessian else h = matrix(NA,3,3)
  out = list(sigma = sigma, skew = skew, kurt = kurt, converge.result = converge.result, hessian = h)
  out
}
