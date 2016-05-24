extract.am.density <-
function(initial.values = rep(NA,10), r, te, s0, market.calls, call.weights = NA, market.puts, put.weights = NA, strikes, lambda = 1, hessian.flag = F, cl = list(maxit=10000))
{

   if ( sum(is.na(call.weights)) >= 1 )  call.weights = rep(1,length(strikes))
   if ( sum(is.na(put.weights))  >= 1 )  put.weights  = rep(1,length(strikes))

  if ( sum(is.na(initial.values)) >= 1 )  
  {
    band      = abs((r - 0.5*(0.3)^2)*te)

    d1 = expand.grid(w.1 = c(0.1, 0.9)  ,w.2 = c(0.1, 0.9),
                     u.1 = c(log(s0) - 5*band, log(s0) - 4 * band) ,u.2 = c(log(s0) - band, log(s0) + band) , u.3 = c(log(s0) + 4* band, log(s0) + 5 * band) , 
                     sigma.1 = c(0.1, 0.4), sigma.2 = c(0.1, 0.4) , sigma.3 = c(0.1, 0.4),
                     p.1 = c(0.1, 0.5) , p.2 = c(0.2, 0.6))

    init.val = apply(d1, 1, mln.am.objective, s0 = s0, r=r, te=te, 
                     market.calls = market.calls, market.puts = market.puts, strikes = strikes, call.weights = call.weights,  put.weights = put.weights, lambda = lambda)

    initial.values = as.numeric(d1[which.min(init.val),])

  }

  optim.obj = optim(initial.values, mln.am.objective, s0 = s0, r=r, te=te, market.calls = market.calls, 
             market.puts = market.puts, strikes = strikes, call.weights = call.weights,  put.weights = put.weights, lambda = 1, hessian = FALSE , control = cl )

  w.1       = optim.obj$par[1] 
  w.2       = optim.obj$par[2]
  u.1       = optim.obj$par[3] 
  u.2       = optim.obj$par[4]
  u.3       = optim.obj$par[5]
  sigma.1   = optim.obj$par[6]
  sigma.2   = optim.obj$par[7]
  sigma.3   = optim.obj$par[8]
  p.1       = optim.obj$par[9]
  p.2       = optim.obj$par[10]

  converge.result = optim.obj$convergence
  if (hessian.flag) h = optim.obj$hessian else h = matrix(NA,10,10)
  out = list(w.1 = w.1, w.2 = w.2, u.1 = u.1, u.2 = u.2, u.3 = u.3, sigma.1 = sigma.1, sigma.2 = sigma.2, sigma.3 = sigma.3, p.1 = p.1, p.2 = p.2, converge.result = converge.result, hessian = h)
  out
}
