#this function is interesting to get the baseline cumulative incidence function ONLY if you have alredy the coefficient estimate
#if not, use the cumulative incidence using the survfit function with a coxph object using the Geskus modified dataset

BBCI<-function( formula , time , status , trans = 1 , cens.code = 0 , data , beta = 0){

  #get the model matrix
  mm<-model.matrix(formula,data)
  
  #get times
  times<-data[,time]

  #modify at risk set
  times2<-times
  times2[!data[,status]%in%c(cens.code,trans)]<-Inf  
  
  #Estimate the survival function of the censoring distribution g_hat
  g1<-with(data,Surv(time=times,event=status==cens.code))
  sg1<-survfit( g1~1 )
  #to avoid errors
  sg1$time[1]<-min(times)

  #compute exp (sum of B'Z )
  beta2<-c(0,beta)
  bz<-t(exp(mm%*%beta2))

  #get tj  (times at which en event of interest occurs)
  nt<-sort(times[data[,status]==trans])

  #get the values of g_hat at tj and at all times   
  ghat<-sapply( nt , function(x) sg1$surv[tail(which( sg1$time <= x),1)])
  ghat2<-sapply( times , function(x) sg1$surv[tail(which( sg1$time <= x),1)])
  
  #get the results in matrix and compute the wheits for the breslow analog for 1-cuminc function
  matg1.1<-matg<-matrix(ghat,nrow=length(times),ncol=length(nt),byrow=T)
  matg2<-matrix(ghat2,nrow=length(times),ncol=length(nt))
  matg3<- matg2 > matg
  matg1.1[matg3]<-matg2[matg3]
  matg4<-matg/matg1.1

  #get the matrix of indicators of times in the modified risk set
  matt4<-matrix(times,nrow=length(times),ncol=length(nt))
  matt<-matrix(times2,nrow=length(times),ncol=length(nt))
  matt2<-matrix(nt,nrow=length(times),ncol=length(nt),byrow=T)
  matt3<-matt >= matt2
  
  #combine indicators en wheits
  matt5<-matt3*matg4
  
  #get the Breslow-like baseline cumulative sub-distribution hazard
  cs<-cumsum(1/(bz%*%matt5))
  
  #get the Breslow-like baseline cumulative incidence function
  s2<- 1 - exp(-cs)  
  
  #get the result in a data frame
  df1<-data.frame(time = nt, est = s2)
  return(df1)
}