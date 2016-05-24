 LK.test<-function(x,type="CO", nsim=1000){
  #Family of the test statistics based on the Laplace 
  LK.statistic<-function(x,type){
      sumKeepAttr <- function( x, keepNames = FALSE, na.rm = FALSE ) {
       xAttr <- attributes( x )
       if( !keepNames ) {
        xAttr$names <- NULL
       }
       x <- sum( x, na.rm = na.rm )
       mostattributes( x ) <- xAttr
       return( x )
    }
  
    TYPE <- deparse(substitute(type))
    n = length(x)
    #EMV
    e = n/sum(x)    
    y=e*x
    ### Test of Cox and Oakes
    if(type=="CO"){
      I=sqrt(n+e*sum(x*(log(y))^2 )-(sum(x*log(y)) )^2 *e^2/n  )
      CO=(n+sum(log(x)-n*x*log(x)/sum(x)))/I
      CO=CO^2
      LK.statistic <- CO
    }### Test of Wald
    else if(type=="W"){
      m <- MLEst(x)
      eta=m$eta
      beta=m$beta
      I=(n/6)*pi^2
      W=I*(beta-1)^2
      LK.statistic <-W
    }### Test of LR
    else if(type=="LR"){  
      m <- MLEst(x)
      eta=m$eta
      beta=m$beta
      #Calcul de la statistique de test
      LR=2*n*log(beta*sum(x)/sum(x^beta))+2*(beta-1)*sum(log(x))
      
      LK.statistic <-LR
    }
    else stop(paste("unknown ", TYPE, "!"))
    return(list(statistic=LK.statistic,lambda=e)) 
  }
  DNAME <- deparse(substitute(x))
  TYPE <- deparse(substitute(type))
  n <- length(x)
  if(sum(x<0)){
    stop(paste("Data ", DNAME, " is not a positive sample"))}
  if(nsim<100){
    warning("small values of Monte-Carlo iterations may affect the value of the p-value")
  }
  if(as.character(type)=="CO"){
    METHOD="Test of Cox and Oakes for the Exponential distribution"
  }  else if(as.character(type)=="W"){
    METHOD="Test of Wald for the Exponential distribution"
  }  else if(as.character(type)=="LR"){
    METHOD="Likelihood ratio test for the Exponential distribution"
  } 
  stat <- LK.statistic(x,type)
  statistic.obs <- stat$statistic
  estimate.obs <- stat$lambda
  fun<-function(y){
    fun <- LK.statistic(y,type)
    return(fun$statistic)
  }
  sim.statistic <- GoFsim(nsim,n,fun) 
  p_val <- sum(sim.statistic>=statistic.obs)/nsim
  LK.test <- list(statistic =c(S=statistic.obs), p.value = p_val, 
                  method=paste(METHOD),
                  estimate = estimate.obs,  
                  data.name=DNAME)
  
  class(LK.test) <-"htest"
  return(LK.test)
}
