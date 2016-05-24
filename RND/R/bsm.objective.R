bsm.objective <-
function(s0, r, te, y, market.calls, call.strikes, call.weights = 1, market.puts, put.strikes, put.weights = 1, lambda = 1, theta)
{
  mu               = theta[1]
  zeta             = theta[2]    
  discount.factor  = exp(-r * te)
  expected.value   = exp(mu + 0.5 * zeta^2)

  ###
  ### Calls
  ###

  d1 = (log(call.strikes) - mu - zeta^2)/zeta
  d2 = (log(call.strikes) - mu)/zeta
  theoretical.calls = discount.factor * ( expected.value * (1 - pnorm(d1))  - call.strikes * (1 - pnorm(d2)) )

  ###
  ### puts
  ###

  c1 = (log(put.strikes) - mu - zeta^2)/zeta
  c2 = (log(put.strikes) - mu)/zeta
  theoretical.puts  = discount.factor * ( put.strikes * pnorm(c2) - pnorm(c1) * expected.value  )
   
  ###
  ### Finally ... the objective function value
  ###

  if ( zeta < 0 ) {obj = 10^7} else { obj = sum(call.weights*(theoretical.calls - market.calls)^2) + sum(put.weights*(theoretical.puts - market.puts)^2) + lambda * (s0*exp(-y*te) - expected.value * discount.factor )^2 }

  obj

}
