transform.par.EXPLN <-
function(this.par,m,fix.mu) {
# Transforms the parameter from its original scale to all real numbers, i.e.
# INPUT:
#    this.par=(p1,p2,...,p_(T-1),mu,sigma,lambda),
# where mu and lambda are typically vectors.   
#
# OUTPUT:
#    this.theta=(w1,w2,...w_(T-1),mu,log(sigma),log(lambda))
# with 
#    w_i = logit((p_1+...+p_i)/(p_1+...+p_(i+1))) for i=1,...,T-1
# and
#    lambda = (lambda_1,lambda_2,...,lambda_m)
#
#
# - this.par is the parameter on the original scale.
# - m is the number of genes analyzed
# - fix.mu is a Boolean indicating whether the parameter mu is kept fixed.
#
# In any case, this.par is full-dimensional, that means, is contains mu independently of 
# the variable fix.mu. The output is lower-dimensional for fix.mu==T.

   
   # determine number of types of cells
   if (length(this.par)==m) {
      TY <- 1
   }
   else {
      TY <- length(this.par)/(m+1)
   }
     
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
   if ((!fix.mu) && (TY>1)) {
      this.mu <- this.par[TY:((m+1)*(TY-1))]
   }
   else {
      this.mu <- NULL
   }
   
   # sigma
   if (TY>1) {
      this.sigma <- this.par[(m+1)*(TY-1)+1]
   }
   else {
      this.sigma <- NULL
   }
   
   # lambda
   if (TY>1) {
      this.lambda <- this.par[(m+1)*(TY-1)+1+(1:m)]
   }
   else {
      this.lambda <- this.par[1:length(this.par)]
   }
   
   # concatenate, transform
   if (TY>1) {
      this.theta <- c(this.w,this.mu,log(this.sigma),log(this.lambda))
   }
   else {
      this.theta <- log(this.lambda)
   }
            
   return(this.theta)
}
