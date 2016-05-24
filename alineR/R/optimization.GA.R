optimization.GA <-
function(Al,data,num,step=5,plot=TRUE){
  para<-random_parameter(num)
  cutoff=num/4
  meanf<-vector()
  fitness<-fitness_count(Al,para=para,Data=data)
  meanf<-mean(fitness)
  for(i in 1:step){
    para<-ES.update(para,fitness,cutoff)
    fitness<-fitness_count(Al=Al,para=para,Data=data)
    meanf<-c(meanf,mean(fitness[1:cutoff]))
    if(plot==TRUE){
      #hist(fitness,main="distribution of fitness")
    }
  }
  if(plot==TRUE){
    plot(meanf,ylim=c(min(meanf)-10,length(Al[1,])),xlab="Number of Interations",ylab="Mean Performance Over All Candidates",main="Convergence For Default Parameter",col="red",)
    abline(h=length(Al[1,]),lty=6,col="blue")
    legend("bottomright",lty=c(6,NA),pch=c(NA,1),col=c("blue","red"),legend=c("Theoretical Optimal","Performnce"))
  }
  z<-list(fitness,para[1:cutoff,])
  names(z)<-c("performance","optimized_parameters")
  return(z)
}
