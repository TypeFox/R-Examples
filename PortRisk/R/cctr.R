cctr <-
function(tickers, weights=rep(1,length(tickers)), start, end, data)
{
  w = weights/sum(weights)
  cctr = w*mctr(tickers, weights, start, end, data)
  return(cctr)
}
