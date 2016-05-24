#####Ploting likelihood process 
LikelihoodHistory=function(Trace,pvalue,Status=FALSE,True.Node) 
{
  rstat=Transfer(pvalue)
  mylog<-0
  like<-0
  for (i in 1:nrow(Trace))
  {
    mu0=mean(rstat[which(Trace[i,]<=0)])
    mu1=mean(rstat[which(Trace[i,]>=1)])
    var0=var(rstat[which(Trace[i,]<=0)])
    var1=var(rstat[which(Trace[i,]>=1)])
    if (is.na(var0)){var0=var1
    }else if (is.na(var1)){var1=var0}
    for(num1 in 1:ncol(Trace)){
      
      if(Trace[i,num1]<=0){
        mylog[num1]<--(rstat[num1]-mu0)^2/(2*var0)-log(sqrt(var0))
      }else {
        mylog[num1]<--(rstat[num1]-mu1)^2/(2*var1)-log(sqrt(var1))
      }
    }
    like=c(like,sum(mylog))
  }
  if (Status==TRUE){
    
    mu0=mean(rstat[which(True.Node<=0)])
    mu1=mean(rstat[which(True.Node>=1)])
    var0=var(rstat[which(True.Node<=0)])
    var1=var(rstat[which(True.Node>=1)])
    if (is.na(var0)){var0=var1
    }else if (is.na(var1)){var1=var0}
    for(num1 in 1:length(rstat)){
      
      if(True.Node[num1]<=0){
        mylog[num1]<--(rstat[num1]-mu0)^2/(2*var0)-log(sqrt(var0))
      }else{
        mylog[num1]<--(rstat[num1]-mu1)^2/(2*var1)-log(sqrt(var1))
      }
    }
    like.true=sum(mylog)
    plot(like[-1],type = "l",ylim=range(like[-1],like.true),xlab="Iteration",ylab="Log-likelihood Value")
    abline(a=like.true,b=0,col="red")
  }else{plot(like[-1],type = "l",xlab="Iteration",ylab="Log-likelihood Value")}
  
  
  invisible(like[-1])
}
