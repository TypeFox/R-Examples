
###EMLR
library("emplik")
library("rootSolve")


### function to calculate Likelihood ratio
# returns numeric value

EMLogLR<-function(x,mean) UseMethod("EMLogLR")
EMLogLR.default<-function(x,mean){
  el.test( x,mean,rep(1,length(mean)),maxit=200, gradtol=1e-7, svdtol = 1e-9, itertrace=FALSE )$'-2LLR'
}
# returns numeric value


##########GAMMA FUNCTION##############
##most simple logliklihood function

exp_factor<-function(x,theta) UseMethod("exp_factor")
exp_factor.default<-function(x,theta){
  con_rate=2
  bart=1
  deltan=1
  # convert x to matrix first
  x=as.matrix(x)
  
  ##con_rate : 1/(L(theta))^(1/(rate*n)), how fast to converge to zero as n goes to inf
  if(length(x[1,])>1) n<-length(x[,1])
  else     n<-length(x)
  #((1/EMLogLR(x,temptheta))^(1/(con_rate*n)))
  # rn<-1/(con_rate*n)
  # gamm<-1+(EMLogLR(x,theta))*rn
  # gamm
  #############SECOND ORDER#############
  rn<-1/(con_rate*n)
  EMLogLR_power<-(bart)*(EMLogLR(x,theta))^(deltan)
  gamm<-1+EMLogLR_power*rn
  gamm                   
}#end exp_factor

# return matrix value


###Qudratic form to solve zeta#######
prime_image<-function(theta_tilda,theta,x) UseMethod("prime_image")
prime_image.default<-function(theta_tilda,theta,x){
  con_rate=2
  bart=1
  deltan=1
  ##length between theta_tilda and theta
  fix_len<- sqrt(sum((theta - theta_tilda)^2))

  ##unit vector of theta_theta_tilda
  if (fix_len!=0){
  unit_vec<- (theta- theta_tilda)/sqrt(sum((theta -  theta_tilda)^2))
  f2<-function(tao){
    result=c()
    for (i in 1: length(tao)){
      result[i]=exp_factor(x,theta_tilda+tao[i]*unit_vec)*tao[i]-fix_len
    }
    return(result)
  }
  tao_prime <- uniroot.all(f2,lower=0,upper=fix_len,maxiter = 200 )
  }else{
    tao_prime=0
  }
 
  
  # choose the value that is closest to theta
  if (length(tao_prime)!=1){
    tao_prime=tao_prime[length(tao_prime)]
  }
  ## compute the theta_prime
  theta_prime<-theta_tilda+tao_prime*unit_vec
  
  theta_prime
}#function done

# return numeric value

EEL<-function(x,theta) UseMethod("EEL")
EEL.default<-function(x,theta){
  con_rate=2
  bart=1
  deltan=1
  stat=c()
  x=as.matrix(x)
  theta_tilda=colMeans(x)
  stat$theta<-theta
  stat$estimating="x-theta"
  stat$prime<-prime_image(theta_tilda=theta_tilda, theta=theta, x)
  
  # the expansion factor is calculated using OEL l(theta)
  
  stat$expansion<-exp_factor(x,theta=stat$prime)
  stat$oel_log<-EMLogLR(x,theta)
  stat$eel_log<-EMLogLR(x,stat$prime)
  stat$call<-match.call()
  
  class(stat)<-"EEL"
  stat
}



print.EEL<-function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\n log eel ratio: \n")
  print(x$eel_log)
}


summary.EEL<-function(object,...)
{
  cat("Call:\n")
  print(object$call)
  cat("\n theta: \n")
  print(object$theta)
  cat("\n estimating equation: \n")
  print(object$estimating)
  cat("\n log oel ratio: \n")
  print(object$oel_log)
  cat("\n prime image: \n")
  print(object$prime)
  cat("\n expasion factor: \n")
  print(object$expansion)
  cat("\n log eel ratio: \n")
  print(object$eel_log)
}



