correct_gamma=function(P_data,rep,n.eta, sample_gamma){
  P=list()
  for (i in 1:rep){
    P[[i]]=P_data[(1+n.eta*(i-1)):(n.eta*i),]
  }
  correct_sample_G=lapply(1:length(sample_gamma), function(i){t(as.matrix(P[[i]]))%*%sample_gamma[[i]]})
  
  correct_gamma=do.call(rbind,lapply(1:length(correct_sample_G), function(i) {cbind(i,correct_sample_G[[i]])}))
  colnames(correct_gamma)[1]=c("replication")
  write.table(correct_gamma, sep=",", file="correct_gamma.csv", row.names=FALSE)
}

