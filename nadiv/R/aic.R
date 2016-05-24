aic <- function(logLik, fp, n = NULL){
   if(is.numeric(n)){
      AICc <- -2*(logLik - fp * (n / (n - fp - 1)))
      delta_AIC <- AICc - min(AICc)
   } else{
        AIC <- -2*(logLik - fp)
        delta_AIC <- AIC - min(AIC)
     } 
   AIClik <- exp(-0.5*delta_AIC)
   w <- AIClik / sum(AIClik)
   
   if(is.numeric(n)){
      return(list(AICc = AICc, delta_AIC = delta_AIC, AIClik = AIClik, w = w))
   } else{
        return(list(AIC = AIC, delta_AIC = delta_AIC, AIClik = AIClik, w = w))
     }
}

