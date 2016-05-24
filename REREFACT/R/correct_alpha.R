correct_alpha=function(P_data,rep, n.eta, sample_alpha){
  P=list()
  for (i in 1:rep){
    P[[i]]=P_data[(1+n.eta*(i-1)):(n.eta*i),]
  }
  
  correct_sample_A=lapply(1:length(sample_alpha), function(i){t(as.matrix(P[[i]]))%*%sample_alpha[[i]]})
  correct_alpha=do.call(rbind,lapply(1:length(correct_sample_A), function(i) {cbind(i,correct_sample_A[[i]])})) 
  colnames(correct_alpha)[1]=c("replication")
  write.table(correct_alpha, sep=",", file=paste("correct_alpha.csv",sep=""), row.names=FALSE)
}


