get.point.estimate <-
function(market.calls, call.strikes, r , te)
{

  len = length(call.strikes) - 1
  tmp.est = numeric(len)

  for (i in 2:len)
  {
    tmp.1      =  2 * exp(r * te) / (call.strikes[i+1] - call.strikes[i-1])
    tmp.2      =  (market.calls[i+1] - market.calls[i])/(call.strikes[i+1] - call.strikes[i])
    tmp.3      =  (market.calls[i-1] - market.calls[i])/(call.strikes[i-1] - call.strikes[i])
    tmp.est[i] =  tmp.1 * ( tmp.2 - tmp.3 )
                                
  }
    point.estimates  = tmp.est[2:len]
    point.estimates  
}
