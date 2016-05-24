AR <-
function(x, params = list(AR = list(p=1))){
 
  p = params$AR$p
  ar.ls <- pacf(x, plot = FALSE)$acf[1:p]
  
  x.sur <- arima.sim(model=list(ar = ar.ls), n = length(x))
  
  return(invisible(x.sur))
}
