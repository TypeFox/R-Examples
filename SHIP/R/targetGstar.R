targetGstar <-
function(x,genegroups) {
  T1   <- target.help(genegroups)
  T2   <- matrix(nrow=length(genegroups),ncol=length(genegroups),data=0)
  corm <- T1*cor(x)
  diag(corm)<-0
  sum.pos <- sum(corm[T1==1 & corm > 0])
  sum.neg <- sum(corm[T1==1 & corm < 0])
   
  cora.pos <- sum.pos/sum(corm>0)
  cora.neg <- sum.neg/sum(corm<0)

  for (i in 1:length(genegroups)) {
    for (j in 1:i) {
      if (i!=j & T1[i,j]==1 & corm[i,j] > 0) T2[i,j] <- cora.pos*(sd(x[,i])*sd(x[,j]))
      if (i!=j & T1[i,j]==1 & corm[i,j] < 0) T2[i,j] <- cora.neg*(sd(x[,i])*sd(x[,j]))
      if (i==j)T2[i,j] <- cov(x[,i],x[,j])
      T2[j,i]<-T2[i,j]
      }
  }
  T2
}

