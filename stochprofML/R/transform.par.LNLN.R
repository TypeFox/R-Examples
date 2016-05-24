transform.par.LNLN <-
function(this.par,m,fix.mu) {
# Transforms the parameter from its original scale to all real numbers, i.e.
# INPUT:
#    this.par=(p1,p2,...,p_(T-1),mu,sigma)   
#
# OUTPUT:
#    this.theta=(w1,w2,...w_(T-1),mu,log(sigma))
# with 
#    w_i = logit((p_1+...+p_i)/(p_1+...+p_(i+1))) for i=1,...,T-1
#
#
# - this.par is the parameter on the original scale.
# - m is the number of genes analyzed
# - fix.mu is a Boolean indicating whether the parameter mu is kept fixed.
#
# In any case, this.par is full-dimensional, that means, is contains mu independently of 
# the variable fix.mu. The output is lower-dimensional for fix.mu==T.

   
   # determine number of types of cells
   TY <- length(this.par)/(m+1)
        
   # (p_1,...,p_(T-1)) --> w
   if (TY>1) {
      this.p <- this.par[1:(TY-1)]
      this.w <- rep(NA,length(this.p))
      if (TY>2) {
         for (i in 1:(length(this.w)-1)) {
            num <- sum(this.p[1:i])
            denom <- sum(this.p[1:(i+1)])
            if (denom!=0) {
               this.w[i] <- stochprof.logit(num/denom)
            }
            else {
               this.w[i] <- stochprof.logit(0)
            }
         }
      }
      this.w[length(this.w)] <- stochprof.logit(sum(this.p))   
   }
   else {
      this.w <- NULL
   }
      
   # mu
   if (!fix.mu) {
      this.mu <- this.par[TY:(length(this.par)-1)]
   }
   else {
      this.mu <- NULL
   }
   
   # sigma
   this.sigma <- this.par[length(this.par)]
   
   # concatenate, transform
   this.theta <- c(this.w,this.mu,log(this.sigma))
            
   return(this.theta)
}
