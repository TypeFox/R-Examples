volatility <-
function(tickers, start, end, data)
{
  p = length(tickers)
  if(p>0)
  {
    dat = access(tickers, start, end, data)
    # Check if sufficient data (at least 2 values needed) is available for all the tickers to compute variance
    if(p==1)
      chk = sum(is.na(dat)) < (length(dat)-1)
    else
      chk = colSums(is.na(dat)) < (nrow(dat)-1)
    # If not, delete columns with insufficient data
    dat = dat[,chk]
    if(sum(chk)!=p)
      warning("insufficient data available for the ticker(s) ", paste(tickers[!chk], collapse=", "), "; at least 2 values for each of the tickers should be available for the given time period", call. = FALSE)
    tickers = tickers[chk]
    p = length(tickers)
  }
  
  if(p==0) stop("function execution stopped due to insufficient data", call. = FALSE)
  
  sigma = rep(NA,p)
  for(i in 1:p)
    sigma[i] = sd(dat[,i], na.rm=TRUE)
  names(sigma) = tickers
  return(100*sigma)
}
