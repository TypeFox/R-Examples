gb.objective <-
function(theta, r, te, y, s0, market.calls, call.strikes, call.weights = 1, market.puts, put.strikes, put.weights = 1, lambda = 1)
{
  a     = theta[1]
  b     = theta[2]   
  v     = theta[3]
  w     = theta[4]  
  discount.factor  = exp(-r * te)
  expected.value   = b * beta(v + 1/a, w - 1/a)/beta(v,w)   ### See page 395 of jpk

  ###
  ### Calls
  ###

  theoretical.calls = price.gb.option(r = r, te = te, s0 = s0, k = call.strikes, y = y, a = a, b = b, v = v, w = w)$call

  ###
  ### Puts
  ###

  theoretical.puts = price.gb.option(r = r, te = te, s0 = s0, k = put.strikes, y = y, a = a, b = b, v = v, w = w)$put
   
  ###
  ### Finally ... the objective function
  ###

  if ( (a <= 0) | (a > 20) | (b <= 0) | (v <= 0) | (w <= 0) | (w <= 4/a) ) { obj = 10^7 } else { obj = sum(call.weights*(theoretical.calls - market.calls)^2) + sum(put.weights*(theoretical.puts - market.puts)^2) + lambda * (s0*exp(-y*te) - expected.value * discount.factor )^2 }
 
  obj

}
