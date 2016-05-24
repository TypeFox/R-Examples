ew.objective <-
function(theta, r, y, te, s0, market.calls, call.strikes,  call.weights = 1, lambda = 1)
{

  sigma       = theta[1]
  skew        = theta[2] 
  kurt        = theta[3]  
      
  discount.factor   =  exp(-r * te)
  m                 =  log(s0) + (r - y - (sigma^2)/2) * te
  expected.value    =  exp(m + 0.5 * sigma^2 * te)
  theoretical.calls =  price.ew.option( r = r, te = te, s0 = s0, k = call.strikes, sigma = sigma, y = y, skew = skew, kurt = kurt)$call

  ###
  ###  Identify the values that may have produced NaN and remove them
  ###

  NaN.index = which(theoretical.calls == "NaN")
  if ( length(NaN.index) >=1 ) {theoretical.calls = theoretical.calls[-NaN.index]; market.calls = market.calls[-NaN.index]}

  ###
  ### Finally ... the objective function value
  ###

  if ( ( sigma < 0 ) | ( kurt < 0) ) {obj = 10^7} else { obj = sum(call.weights * (theoretical.calls - market.calls)^2) + lambda * (s0 * exp(-y*te) - expected.value * discount.factor )^2 }

  obj

}
