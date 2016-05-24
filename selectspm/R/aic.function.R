aic.function <-
function(r,dtheta, npar){
      
      n<-length(r)
       # only to "see" the formulas as in   Burhnm y Anderson 2002,
       # we rename dtheta as RSS  and npar as K
      RSS <- dtheta
      K <- npar
     
      #loglikelihood
      LL <- (-n/2)*log(RSS/n)
      AIC <- (-2*LL )+ (2*K)
      AICc <- AIC +( (2*K*(K+1))/(n-K-1) )
      return(data.frame(n=n, K=K, RSS=RSS, LL=LL, AIC=AIC,AICc=AICc))
   }
