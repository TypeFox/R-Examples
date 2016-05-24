atkinson.GB2 <-
function(b,a,p,q,epsilon=NULL,ylim=c(0,1000000),zeroapprox=0.01){
  if(is.null(epsilon)) epsilon<-1
  alpha<- 1 - epsilon # see Cowell (2000) p.115
  GE<-entropy.GB2(b,a,p,q,alpha,ylim,zeroapprox)
  if(epsilon<0){
    print("epsilon must be greater or equal to zero")
  }else if(epsilon==1){
    A<-1-exp(-GE)
  }else{
    A<-1-(1+epsilon*(epsilon-1)*GE)^(1/(1-epsilon)) #from A note on the relationship between Atkinson index... , note: 1st transformation wrong!  
  }
    
  return(A)
}
