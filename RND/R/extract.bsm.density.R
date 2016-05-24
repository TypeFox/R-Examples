extract.bsm.density <-
function(initial.values = c(NA, NA), r, y, te, s0, market.calls, call.strikes, call.weights = 1, market.puts, put.strikes, put.weights = 1, lambda = 1, hessian.flag = F, cl = list(maxit=10000))
{
  if ( sum(is.na(initial.values)) >= 1 )  
  {
    band     = (r - y - 0.5*(0.3)^2)*te
    bsm.grid = expand.grid(mu = seq(from = log(s0) - band, to = log(s0) + band, length.out = 10), zeta = sqrt(te) * seq(from = 0.05, to = 0.9, length.out = 10))
    bsm.vals = apply(bsm.grid, 1, bsm.objective, r = r, y = y, te = te, s0 = s0, market.calls = market.calls, call.strikes = call.strikes, 
                             market.puts = market.puts, put.strikes = put.strikes, lambda = 1, call.weights = call.weights,  put.weights = put.weights)
    initial.values = as.numeric(bsm.grid[which.min(bsm.vals),])
  }

  optim.obj = optim(initial.values, bsm.objective, r=r, y = y, te=te, s0 = s0, market.calls = market.calls, 
               call.strikes = call.strikes, market.puts = market.puts, put.strikes = put.strikes, lambda = lambda, hessian = hessian.flag , 
               call.weights = call.weights,  put.weights = put.weights, control=cl )

  mu   = optim.obj$par[1]
  zeta = optim.obj$par[2]
  converge.result = optim.obj$convergence
  if (hessian.flag) h = optim.obj$hessian else h = matrix(NA,2,2)
  out = list(mu = mu, zeta = zeta, converge.result = converge.result, hessian = h)
  out
}
