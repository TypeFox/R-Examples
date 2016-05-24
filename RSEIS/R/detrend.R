`detrend` <-
function(x)
  {
    ##  remove linear trend from a time series
   
    n = length(x)
    
    t = seq(1,n)

    r = lm(x ~ t)
    b = r$coefficients[1]
    m = r$coefficients[2]
    
    newf = x - (m*t+b)
    
    return(newf)
}

