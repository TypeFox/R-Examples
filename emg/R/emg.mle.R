emg.mle <- function(x, lower=NA, upper=NA)
{
   # To provide some reasonable defaults
   if(is.na(lower))
   {
     lower <- list(mu=min(x), sigma=sd(x)/10, lambda=0.001/mean(x))
   }
   if(is.na(upper))
   {
     upper <- list(mu=max(x), sigma=(max(x)-min(x))/4, lambda=100/mean(x))
   }
  
   mle(function(mu, sigma, lambda){
     #print(paste(x, mu, sigma, lambda))
     emg.nllik(x, mu, sigma, lambda)},
             method='L-BFGS-B',
             lower=lower,
             upper=upper,
             start=list(mu=median(x), sigma=sd(x), lambda=1/mean(x)))
}