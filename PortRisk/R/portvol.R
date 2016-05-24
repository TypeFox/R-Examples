portvol <-
function(tickers, weights = rep(1,length(tickers)), start, end, data)
{
  p = length(tickers)
  if(length(weights)!=p)
    stop("unequal lengths of the arguments 'tickers' & 'weights'", call. = FALSE)
  
  if(p>1)
  {
    dat = access(tickers, start, end, data)
    # Check if sufficient data (at least 2 values needed) is available for all the tickers to compute variance-covariance matrix
    chk = colSums(is.na(dat)) < (nrow(dat)-1)
    # If not, delete columns with insufficient data
    dat = dat[,chk]
    if(sum(chk)!=p)
      warning("insufficient data available for the ticker(s) ", paste(tickers[!chk], collapse=", "), "; at least 2 values for each of the tickers should be available for the given time period", call. = FALSE)
    tickers = tickers[chk]
    p = length(tickers)
  }
  
  if(p<=1) return(volatility(tickers, start, end, data))
  
  Sigma = cov(dat, use="pairwise.complete.obs")
  w = weights/sum(weights) # weights
  sigma_p = (t(w) %*% Sigma %*% w)[1,1]
  return(100*sqrt(sigma_p))
}
