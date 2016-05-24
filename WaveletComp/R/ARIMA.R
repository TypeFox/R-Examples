ARIMA <-
function(x, params=list(ARIMA=list(p=1, q=1, include.mean=T, sd.fac=1, trim = F, trim.prop = 0.01 ))) {

   n = length(x)
   
   p = params$ARIMA$p                        # ar order
   q = params$ARIMA$q                        # ma order
   include.mean = params$ARIMA$include.mean  # include intercept?
   
   sd.fac = params$ARIMA$sd.fac              # magnification factor to the residual standard deviation in simulation 
   trim = params$ARIMA$trim                  # simulate trimmed data?
   trim.prop = params$ARIMA$trim.prop        # high/low trimming proportion
      
   my.arima = arima(x, order = c(p,0,q), include.mean=include.mean)
     
   ar.coef = 0
   ma.coef = 0
   intercept = 0
     
   if (p>0) { 
       ar.coef = as.vector(my.arima$coef[1:p])
   }    
   if (q>0) {
       ma.coef = as.vector(my.arima$coef[(p+1):(p+q)])         
   }
   
   if (include.mean==T) {
       intercept = as.numeric(my.arima$coef[p+q+1])
   }
   
   sd = sqrt(my.arima$sigma2)
   
   
   if (trim == F) {
   
       eps = rnorm(n, sd=sd.fac*sd)
       x.sur <- arima.sim(model=list(ar=ar.coef, ma=ma.coef), n = n, innov=eps) + intercept
       
   }
   
   if (trim == T) {
   
       q01 = quantile(x, probs=trim.prop) - intercept
       q99 = quantile(x, probs=1-trim.prop) - intercept
       
       eps = rnorm(n+p+q, sd=sd.fac*sd)
       x.sur = rep(0,max(p,q))
       
       for (i in ((max(p,q)+1):(n+p+q))) {
            x.sur[i] = eps[i]
            if (p>0) {
                x.sur[i] = x.sur[i] + sum(ar.coef*x.sur[(i-1):(i-p)])
            }
            if (q>0) {
                x.sur[i] = x.sur[i] + sum(ma.coef*eps[(i-1):(i-q)])
            }
            if (x.sur[i] < q01) { x.sur[i] = q01 } 
            if (x.sur[i] > q99) { x.sur[i] = q99 }
        }
        x.sur = x.sur + intercept
        x.sur = x.sur[(p+q+1):(n+p+q)]
        
   }
   
   return(invisible(x.sur))
}
