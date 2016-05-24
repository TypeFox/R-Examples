#Source: rtn.R
#random variates from truncated normal distribution
rtn<-function(n, zmu=0, zsig=1, c=1.96, side = c("left", "right")){
  side <- match.arg(side)
  u <- runif(n)
  if(side=="left") 
    y <- -qnorm(u*pnorm(-c, mean=-zmu, sd=zsig),mean=-zmu,sd=zsig)
  else 
    y<-qnorm(u*pnorm(c, mean=zmu, sd=zsig),mean=zmu,sd=zsig)
  y
} 

