correct_psi=function(P_data,rep,n.eta, sample_psi){
  P=list()
  for (i in 1:rep){
    P[[i]]=P_data[(1+n.eta*(i-1)):(n.eta*i),]
  }
  
  correct_sample_psi=lapply(1:length(sample_psi), function(i){t(as.matrix(P[[i]]))%*%sample_psi[[i]]%*%as.matrix(P[[i]])})
  
  correct_psi=do.call(rbind,lapply(1:length(correct_sample_psi), function(i) {cbind(i,correct_sample_psi[[i]])}))
  colnames(correct_psi)[1]=c("replication")
  write.table(correct_psi, sep=",", file="correct_psi.csv", row.names=FALSE)
 
}


