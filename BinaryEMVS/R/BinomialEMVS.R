#' Variable Selection For Binary Data Using The EM Algorithm
#'
#' Conducts EMVS analysis
#'
#' @param y responses in 0-1 coding
#' @param x X matrix
#' @param type probit or logit model
#' @param epsilon tuning parameter
#' @param v0s tuning parameter, can be vector
#' @param nu.1 tuning parameter
#' @param nu.gam tuning parameter
#' @param lambda.var tuning parameter
#' @param a tuning parameter
#' @param b tuning parameter
#' @param beta.initial starting values
#' @param sigma.initial starting value
#' @param theta.inital startng value
#' @param temp not sure
#' @param p not sure
#' @param n not sure
#' @param SDCD.length not sure
#' 
#' @return probs is posterior probabilities
#'
#' @examples
#' #Generate data
#' set.seed(1)
#' n=25;p=500;pr=10;cor=.6
#' X=data.sim(n,p,pr,cor)
#' 
#' #Randomly generate related beta coefficnets from U(-1,1)
#' beta.Vec=rep(0,times=p)
#' beta.Vec[1:pr]=runif(pr,-1,1)
#' 
#' y=scale(X%*%beta.Vec+rnorm(n,0,sd=sqrt(3)),center=TRUE,scale=FALSE)
#' prob=1/(1+exp(-y))
#' y.bin=t(t(ifelse(rbinom(n,1,prob)>0,1,0)))
#' 
#' result.probit=BinomialEMVS(y=y.bin,x=X,type="probit")
#' result.logit=BinomialEMVS(y=y.bin,x=X,type="logit")
#' 
#' which(result.probit$posts>.5)
#' which(result.logit$posts>.5)
#'
#' @export
BinomialEMVS=function(y,x,type="probit",epsilon=.0005,
                       v0s=ifelse(type=="probit",.025,5),
                       nu.1=ifelse(type=="probit",100,1000),
                       nu.gam=1,lambda.var=.001,a=1,b=ncol(x),
                       beta.initial=NULL,
                       sigma.initial=1,theta.inital=.5,temp=1,p=ncol(x),n=nrow(x),SDCD.length=50){
  
  if(type=="probit")
  {
    result=
      EMVS.probit(y=y,x=x,epsilon=epsilon,v0s=v0s,nu.1=nu.1,nu.gam=nu.gam,a=a,b=b,
                  beta.initial=beta.initial,sigma.initial=sigma.initial,theta.inital=.5,
                  temp=temp,p=p,n=n)
  }
  
  if(type=="logit")
  {
    y[y==0]=-1
    result=
      EMVS.logit(y=y,x=x,epsilon=epsilon,v0s=v0s,nu.1=nu.1,nu.gam=nu.gam,a=a,b=b,
                 beta.initial=beta.initial,sigma.initial=sigma.initial,theta.inital=.5,
                 temp=temp,p=p,n=n,lambda.var=lambda.var,SDCD.length=SDCD.length)
  }
  return(result)  
}