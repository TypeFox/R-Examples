targetCor <-
function(x,genegroups) {
         
  T1 <- target.help(genegroups)
  T2 <- matrix(nrow=length(genegroups),ncol=length(genegroups),data=0)     

  for (i in 2:length(genegroups)) {
     for (j in 1:(i-1)) if (T1[i,j]==1) T1[i,j] <- ifelse(cor.test(x[,i],x[,j])$p.value < 0.05,1,0)
  }
  
  corm <- T1*cor(x)     
  diag(corm) <- 0  
  
  cora <- ifelse(sum(corm!=0)==0,0,sum(colSums(corm))/sum(corm!=0))
  
  for (i in 1:length(genegroups)) {
    for (j in 1:i) {
      if (i!=j & T1[i,j]==1) T2[i,j] <- cora*(sd(x[,i])*sd(x[,j]))
      if (i==j)T2[i,j] <- cov(x[,i],x[,j])
      T2[j,i]<-T2[i,j]
      }
  }
  T2
}

