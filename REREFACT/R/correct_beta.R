correct_beta=function(P_data, rep, n.eta, sample_beta){
  P=list()
  for (i in 1:rep){
    P[[i]]=P_data[(1+n.eta*(i-1)):(n.eta*i),]
  }
  correct_sample_B=lapply(1:length(sample_beta), function(i){t(as.matrix(P[[i]]))%*%sample_beta[[i]]%*%as.matrix(P[[i]])})
  
  correct_beta=do.call(rbind,lapply(1:length(correct_sample_B), function(i) {cbind(i,correct_sample_B[[i]])}))
  colnames(correct_beta)[1]=c("replication")
  write.table(correct_beta, sep=",", file="correct_beta.csv", row.names=FALSE)
}

