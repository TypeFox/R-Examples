extract.mln.density <-
function(initial.values = c(NA,NA,NA,NA,NA), r, y, te, s0, market.calls, call.strikes, call.weights = 1, market.puts, put.strikes, put.weights = 1, lambda = 1, hessian.flag = F, cl = list(maxit=10000))
{

  if ( sum(is.na(initial.values)) >= 1 )  
  {
    band      = (r - y - 0.5*(0.3)^2)*te
    mln.grid  = expand.grid(alpha.1 = seq(from = 0.1, to = 0.90, by = 0.05), meanlog.1 = seq(from = log(s0) - band, to = log(s0) + band, length.out = 4),
                                 meanlog.2 = seq(from = log(s0) - band, to = log(s0) + band, length.out = 4), sdlog.1 = sqrt(te) * seq(from = 0.05, to = 0.9, length.out = 7),
                                 sdlog.2 = sqrt(te) * seq(from = 0.05, to = 0.9, length.out = 7))
    mln.vals  = apply(mln.grid, 1, mln.objective, r = r, y = y, te = te, s0 = s0, market.calls = market.calls, call.strikes = call.strikes, 
                             market.puts = market.puts, put.strikes = put.strikes, call.weights = call.weights,  put.weights = put.weights, lambda = 1)
    initial.values = as.numeric(mln.grid[which.min(mln.vals),])
  }

  optim.obj = optim(initial.values, mln.objective, r=r, y=y, te=te, s0 = s0, market.calls = market.calls, 
               call.strikes = call.strikes, market.puts = market.puts, put.strikes = put.strikes, lambda = lambda, hessian = hessian.flag , 
               call.weights = call.weights,  put.weights = put.weights, control = cl )

  alpha.1     = optim.obj$par[1]
  meanlog.1   = optim.obj$par[2]
  meanlog.2   = optim.obj$par[3]
  sdlog.1     = optim.obj$par[4]
  sdlog.2     = optim.obj$par[5]

  converge.result = optim.obj$convergence
  if (hessian.flag) h = optim.obj$hessian else h = matrix(NA,5,5)
  out = list(alpha.1 = alpha.1, meanlog.1 = meanlog.1, meanlog.2 = meanlog.2, sdlog.1 = sdlog.1, sdlog.2 = sdlog.2, converge.result = converge.result, hessian = h)
  out
}
