extract.shimko.density <-
function(market.calls, call.strikes, r, y, te, s0, lower, upper)
{

  implied.volatilities = compute.implied.volatility(r = r, te = te, s0 = s0, k = call.strikes, 
                                                    y = y, call.price = market.calls, lower = lower , upper = upper)

  implied.curve.obj = fit.implied.volatility.curve(x = implied.volatilities,k = call.strikes)


  shimko.density = dshimko(r = r, te = te, s0 = s0, k = call.strikes, y = y, 
                           a0 = implied.curve.obj$a0, a1 = implied.curve.obj$a1, a2 = implied.curve.obj$a2)

  out = list(implied.curve.obj = implied.curve.obj, shimko.density = shimko.density, implied.volatilities = implied.volatilities)
  out

}
