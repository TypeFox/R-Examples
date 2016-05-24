#####Iteration3_Mclust
Iteration3_Mclust<-function(iter,wholeindex,hodcmclust,net,pirhopair,choice,rstat,show.steps,showlikelihood,likelihood.frequency){
  z<-wholeindex
  total<-matrix(rep(0,length(wholeindex)*iter),ncol=length(wholeindex),nrow=iter)
  
  for(jj in 1: iter){
    
    if(jj%%show.steps==0){
      cat("iter: ",jj,"\n")
      flush.console()
    }
    for(num1 in 1:length(wholeindex)){
      ztemp=c()
      pro1<-0
      pro0<-0
      idx = which(net[num1,]==1)
      pro0=sum((z[idx]==0))
      pro1=sum((z[idx]==1))
      if (length(hodcmclust$variance)==1){hodcmclust$variance=c(hodcmclust$variance,hodcmclust$variance)}
      
      log0<-log(pirhopair$pi0[choice])+2*pirhopair$rho0[choice]*pro0-(rstat[num1]-hodcmclust$mean[1])^2/(2*hodcmclust$variance[1])-log(sqrt(2*pi*hodcmclust$variance[1]))
      log1<-log(1-pirhopair$pi0[choice])+2*pirhopair$rho1[choice]*pro1-(rstat[num1]-hodcmclust$mean[2])^2/(2*hodcmclust$variance[2])-log(sqrt(2*pi*hodcmclust$variance[2]))
      
      
      p0<-1/(1+exp(log1-log0))
      p1<-1-p0
      
      if (is.na(p0) | is.na(p1)) {z[num1]=sample(c(1,0),1,prob=c(0.5,0.5))###in case the mu is null
      }else{
        z[num1]<-sample(c(0,1),1,prob=c(p0,p1))}
      
      
      total[jj,num1]<-z[num1]
      
      #ratio[jj,num1]<-total[jj,num1]/jj
      
    }
    
    if(showlikelihood==TRUE){
      if(jj%%likelihood.frequency==0){
        mylog<-0
        mu0=mean(rstat[which(total[jj,]<=0)])
        mu1=mean(rstat[which(total[jj,]>=1)])
        var0=var(rstat[which(total[jj,]<=0)])
        var1=var(rstat[which(total[jj,]>=1)])
        for(num1 in 1:length(rstat)){
          
          if(total[jj,num1]<=0){
            mylog[num1]<--(rstat[num1]-mu0)^2/(2*var0)-log(sqrt(var0))
          }else{
            mylog[num1]<--(rstat[num1]-mu1)^2/(2*var1)-log(sqrt(var1))
          }
          
        }
        cat("Now for the step:" ,jj, "the log-likelihood value is" ,sum(mylog) , "\n")
        flush.console()
      }
    }
  }
  
  return(total)
}
