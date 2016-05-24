correct_lambda=function(P_data,rep,n.eta, sample_lambda){
  P=list()
  for (i in 1:rep){
    P[[i]]=P_data[(1+n.eta*(i-1)):(n.eta*i),]
  }
  correct_sample_lambda=lapply(1:length(sample_lambda), function(i) {sample_lambda[[i]]%*%as.matrix(P[[i]])})
  correct_lambda=do.call(rbind,lapply(1:length(correct_sample_lambda), function(i) {cbind(i,correct_sample_lambda[[i]])}))
  colnames(correct_lambda)[1]=c("replication")
  write.table(correct_lambda, sep=",", file="correct_lambda.csv", row.names=FALSE)
}



