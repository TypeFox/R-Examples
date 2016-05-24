ComputeLambdaAndNuHat.shiftCMPest <-
function(y,lambdainit=mean(y)-min(y),nuinit=1,a=min(y),max=100){
  
  minusloglike <- function(par){
    
    -(sum(y-a)*log(par[1]) - par[2]*sum(log(factorial(y-a))) - length(y)*log(computez.lambdaest(par[1],par[2],max)))}
  
  LambdaNuEst <- nlminb(start=c(lambdainit,nuinit),minusloglike,lower = c(0,0), upper = c(Inf,Inf))
  
  return(LambdaNuEst)
}
