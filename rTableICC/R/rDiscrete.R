rDiscrete <-
function(n=1,pf){
    if ((is.finite(pf)==FALSE) | (is.null(pf)==TRUE) | (min(pf)<0) | (max(pf)>1) | (round(sum(pf))!=1)){
       stop("Sum of probabilities is not equal to one! Please enter a probability function.")
    }
    if ((is.finite(n)==FALSE) | (is.null(n)==TRUE) | (min(n)<0) | (length(n)!=1) ){
      stop("Number of observations must be entered as a finite and positive scalar!")
    }
    
    cdf.initial=t(cumsum(pf))
    r=0
    r=runif(n,0,1)
    cdf=cbind(0,cdf.initial)
    names=array("",(length(cdf.initial)+1))
    rDiscrete=0
    for (i in 1:n){
      for (j in 2:(length(cdf)+1)){
        if ((r[i]>cdf[j-1]) & (r[i]<=cdf[j])){
          rDiscrete[i]=j-1
        }
      }
    }
    cdf.initial=as.array(c(0,cdf.initial))
    list(rDiscrete=rDiscrete,cdf=cdf.initial)
}
