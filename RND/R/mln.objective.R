mln.objective <-
function(theta, r, y, te, s0, market.calls, call.strikes, call.weights=1, market.puts, put.strikes, put.weights=1, lambda = 1)
{
  alpha.1     = theta[1]
  meanlog.1   = theta[2] 
  meanlog.2   = theta[3]
  sdlog.1     = theta[4]
  sdlog.2     = theta[5]  

  discount.factor = exp(-r * te)
  alpha.2         = 1 - alpha.1
  expected.value  = alpha.1 * exp(meanlog.1 + (0.5)*(sdlog.1^2)) + alpha.2 * exp(meanlog.2 + (0.5)*(sdlog.2^2))

  ###
  ### Calls
  ###

  theoretical.calls = price.mln.option(r = r, te = te, y = y, k = call.strikes, alpha.1 = alpha.1, 
                                       meanlog.1 = meanlog.1, meanlog.2 = meanlog.2, sdlog.1 = sdlog.1, sdlog.2 = sdlog.2)$call

  ###
  ### puts
  ###

  theoretical.puts = price.mln.option(r = r, te = te, y = y, k = put.strikes, alpha.1 = alpha.1, 
                                       meanlog.1 = meanlog.1, meanlog.2 = meanlog.2, sdlog.1 = sdlog.1, sdlog.2 = sdlog.2)$put


  ###
  ### Finally ... the objective function
  ###

  if ( (alpha.1 < 0) | (alpha.1 > 1) | (sdlog.1 < 0) | (sdlog.2 < 0) ) {obj = 10^7} else { obj = sum(call.weights * (theoretical.calls - market.calls)^2) + sum(put.weights*(theoretical.puts - market.puts)^2) + lambda * (s0 - exp(y*te) * expected.value * discount.factor )^2 }
  obj

}
