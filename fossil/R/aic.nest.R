`aic.nest` <-
function(comm1, comm2, base=exp(1)) {
  if (length(comm1)!=length(comm2)) {
    warning("sample vectors cannot be different lengths")
    return()
  }
  sum1<-sum(comm1)
  sum2<-sum(comm2)
  p1<-comm1/sum1
  p2<-comm2/sum2
  pc<-(sum1*p1+sum2*p2)/(sum1+sum2)
  comb<-comm1+comm2
  comb<-comb[comb>0]
  Sc<-length(comb[comb>0])
  shan<- -(sum(pc*log(pc[pc>0], base)))
  Lc<-(sum(comm1)+sum(comm2))*shan
  shan1<- -(sum(p1[p1>0]*log(p1[p1>0], base)))
  shan2<- -(sum(p2[p2>0]*log(p2[p2>0], base)))
  L2<-sum(comm1)*shan1 + sum(comm2)*shan2
  aicc<-2*(Lc+Sc)
  aic2<-2*(L2+2*Sc)
  results<-c(aicc, aic2)
  names(results)<-c('AIC.same', 'AIC.different')
  return(results)
}



