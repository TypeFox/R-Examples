extract.gb.density <-
function(initial.values = c(NA,NA,NA,NA), r, te, y, s0, market.calls, call.strikes, call.weights = 1, market.puts, put.strikes, put.weights = 1, lambda = 1, hessian.flag = F, cl = list(maxit=10000))
{

  if ( sum(is.na(initial.values)) >= 1 )  
  {
    gb.grid  = expand.grid(a = seq(from = 0, to = 10, length.out = 10), b = seq(from = s0 - 5, to = s0 + 5,length.out = 10), 
                                  v = seq(from = 0.1, to = 5, length.out = 10), w = seq(from = 0.1, to = 5, length.out = 10))
    bad.ones = gb.grid$a * gb.grid$w
    gb.grid  = gb.grid[-which(bad.ones <= 4),]
    gb.vals  = apply(gb.grid, 1, gb.objective, r = r, te = te, y = y, s0 = s0, market.calls = market.calls, call.strikes = call.strikes, 
                            market.puts = market.puts, put.strikes = put.strikes, call.weights = call.weights,  put.weights = put.weights, lambda = 1)
    initial.values = as.numeric(gb.grid[which.min(gb.vals),])
  }

  optim.obj = optim(initial.values, gb.objective, r=r, te=te, y = y, s0 = s0, market.calls = market.calls, 
               call.strikes = call.strikes, market.puts = market.puts, put.strikes = put.strikes, lambda = lambda, hessian = hessian.flag , 
               call.weights = call.weights,  put.weights = put.weights, control = cl )

  a   = optim.obj$par[1]
  b   = optim.obj$par[2]
  v   = optim.obj$par[3]
  w   = optim.obj$par[4]
  converge.result = optim.obj$convergence
  if (hessian.flag) h = optim.obj$hessian else h = matrix(NA,4,4)
  out = list(a = a, b = b, v = v, w = w, converge.result = converge.result, hessian = h)
  out
}
