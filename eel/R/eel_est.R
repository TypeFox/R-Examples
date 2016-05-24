##########GAMMA FUNCTION##############
##most simple logliklihood function

exp_factor_est<-function(x,theta,equation) UseMethod("exp_factor_est")
exp_factor_est.default<-function(x,theta,equation){
  con_rate=2
  bart=1
  deltan=1
  # convert x to matrix first
  x=as.matrix(x)
  theta=as.vector(theta)
  if(length(x[1,])>1) n<-length(x[,1])
  else     n<-length(x)
  #((1/EMLogLR(x,temptheta))^(1/(con_rate*n)))
  # rn<-1/(con_rate*n)
  # gamm<-1+(EMLogLR(x,theta))*rn
  # gamm
  #############SECOND ORDER#############
  rn<-1/(con_rate*n)
  gx=eval(parse(text=equation))
  gx=as.matrix(gx)
  EMLogLR_power<-(bart)*(EMLogLR(gx,c(rep(0,ncol(gx)))))^(deltan)
  gamm<-1+EMLogLR_power*rn
  gamm                   
}#end exp_factor

# return matrix value


###Qudratic form to solve zeta#######
prime_image_est<-function(theta_tilda,theta,x,equation) UseMethod("prime_image_est")
prime_image_est.default<-function(theta_tilda,theta,x,equation){
  con_rate=2
  bart=1
  deltan=1
  ##length between theta_tilda and theta
  fix_len<- sqrt(sum((theta - theta_tilda)^2))
  
  if (fix_len==0) {
    tao_prime=0
  }
  else{

    unit_vec<- (theta- theta_tilda)/sqrt(sum((theta -  theta_tilda)^2))
  
  ## function of h(zeta)=gmma(n,l(zeta))*zeta
  
  f2<-function(tao){
    result=c()
    for (i in 1: length(tao)){
      result[i]=exp_factor_est(x,theta_tilda+tao[i]*unit_vec,equation)*tao[i]-fix_len
    }
    return(result)
  }
  ## solve the zeta 
 
    tao_prime <- uniroot.all(f2,lower=0,upper=fix_len,maxiter = 200 )
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



EEL_est<-function(x,theta, theta_tilda, equation) UseMethod("EEL_est")
EEL_est.default<-function(x,theta, theta_tilda,equation){
  x=as.matrix(x)
  theta=as.matrix(theta)
  con_rate=2
  bart=1
  deltan=1
  
  stat=c()
  stat$theta<-theta
  stat$prime<-prime_image_est(theta_tilda=theta_tilda, theta=theta, x,equation)
  # the expansion factor is calculated using OEL l(theta)
  stat$estimating<-equation
  stat$expansion<-exp_factor_est(x,theta=stat$prime,equation)
  stat$oel_log<-EMLogLR(x,theta)
  stat$eel_log<-EMLogLR(x,stat$prime)
  stat$call<-match.call()
  
  class(stat)<-"EEL"
  stat
}







