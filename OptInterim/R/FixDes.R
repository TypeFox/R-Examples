"FixDes" <-
function(B.init,m.init,alpha,beta,param,x,num.arm,r=0.5)
{  
  shape0 <- param[1]
  scale0 <- param[2]
  shape1 <- param[3]
  scale1 <- param[4]

  if(!identical(num.arm,1) && !identical(num.arm,2))
  {
    stop("the number of treatment arms must be either one or two")
  }   
 
  if(identical(num.arm,1) && !identical(r,0.5))
  {
    stop("randomization ratio not equal to 0.5 for a single-arm study")
  }
  
  if(identical(num.arm,2) && r*(1-r)<=0)
  {
    stop("invalid randomization ratio r")
  }
    
  if(!identical(length(B.init),length(m.init)))
  {
    stop("projected patient times and numbers should be of equal length")
  }

  if(x>=max(B.init))
  {
    stop("the survival time of interest is beyond the projected accrual time")
  }
  
  if(qnorm(1-alpha)+qnorm(1-beta)<=0)stop("Stopped because alpha, beta are invalid")

  p0<-s(shape0,scale0,x)
  p1<-s(shape1,scale1,x)
  lam0<-lambda(shape0,scale0,x)
  lam1<-lambda(shape1,scale1,x)

  if(p0>=p1)stop("The null event-free rate exceeds the alternative rate")


  ### sample size and study times based on normal approximation to log
  ### hazard function
  sig21 <- sqrt((1-p1)/p1)  ### delta method applied to binomial variance (based on Taylor expansion) to estimate sig21 in equation 9
  sig20 <- sqrt((1-p0)/p0)
  
  if(identical(num.arm,1)) {
      n0 <-(sig21*(qnorm(1-alpha) + qnorm(1-beta))/((log(lam0)-
            log(lam1))*lam1))^2
      n0 <- ceiling(n0)
      if(n0>floor(sum(m.init)))stop("Sample size exceeds specified accrual rates/time")

      da<-compMDA(B.init,m.init,n0)[1]
  }
  
  if(identical(num.arm,2)) {
      v0 <- (sig20^2)/((1-r)*lam0^2)
      v1 <- (sig21^2)/(r*lam1^2)
      n0 <- (sqrt(v0+v1)*(qnorm(1-alpha) + qnorm(1-beta))/((log(lam0)-
            log(lam1))))^2
      n0 <- ceiling(n0)
      if(n0>floor(sum(m.init)))stop("Sample size exceeds specified accrual rates/time")

      da<-compMDA(B.init,m.init,n0)[1]
  }
  
  ### sample size and study times based on exact calculations (one-arm: exact binomial; two-arm: Fisher exact)

  if(identical(num.arm,1)) {
    n0E<-single.exact(n0,alpha,beta,p0,p1)
    daE<-compMDA(B.init,m.init,n0E)[1]
  }
  
  
  
  if(identical(num.arm,2)) {
    ## package clinfun Fisher Exact test
    n0E<-sum(fe.ssize(p0,p1,alpha*2,1-beta,r/(1-r))["Fisher Exact",1:2])
    daE<-compMDA(B.init,m.init,n0E)[1]
  }  

  

  return(list(n0=n0,DA=da,SL=da+x,n0E=n0E,DAE=daE,SLE=daE+x,C=qnorm(1-alpha)))
}  


##########################################################################################
# for one-arm:                                                                          ##
#    normal approximation is more conservative than exact binomial (n.norm>n.binom)     ##
# for two-arm:                                                                          ##
#    fisher exact test is more conservative than normal approximation (n.fisher>n.norm) ##
#                                                                                       ##
# Fisher exact test is a very conservative test when margins are not fixed              ##
##########################################################################################

