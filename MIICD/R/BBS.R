

BBS<-function( formula , time , status , cens.code = 0 , data , beta = 0 ){

  #get the model matrix
  mm<-model.matrix(formula,data)
  
  #get times
  times<-data[,time]

  
  #compute exp (sum of B'Z )
  beta2<-c(0,beta)
  bz<-t(exp(mm%*%beta2))

  #get tj  (times at which en event of interest occurs)
  nt<-sort(times[data[,status]!=cens.code])
  
  #get the matrix of indicators of times in the modified risk set
  matt4<-matrix(times,nrow=length(times),ncol=length(nt))
  matt2<-matrix(nt,nrow=length(times),ncol=length(nt),byrow=T)
  matt3<-matt4 >= matt2
  
  #get the Breslow-like baseline cumulative sub-distribution hazard
  cs<-cumsum(1/(bz%*%matt3))
  
  #get the Breslow-like baseline cumulative incidence function
  s2<- exp(-cs)  
  
  #get the result in a data frame
  df1<-data.frame(time = nt, surv = s2)
  return(df1)
}