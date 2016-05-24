rerefact=function(n.factor, n.eta, n.var, pop_lambda, sample_lambda) 
{ 
  set.seed(1234) # set.seed is used to control for the order (index) of possible permutation. 
  s_n.factor=permn(n.factor) #s_n.factor is the total number of permutation by the order of n.factor. # s_n.factor
  perm_order=length(s_n.factor) # perm_order
  v_sign=as.matrix(expand.grid(rep(list(c(1,-1)), n.factor))) # v_sign is the possible of n.factor permuting sign. # v_sign
  colnames(v_sign)=paste("sign", 1:n.factor, sep="")
  
  ###############################################################################################
  # Step #1 
  ###############################################################################################
  
  ###############################################################################################
  #1. Step 1 determines the total number of equivalent forms, I, of the vector of factors, eta: n.perm.
  ###############################################################################################
  
  n.perm=length(s_n.factor)*nrow(v_sign) # Total number of possible permutations: n.perm
  
  ###############################################################################################
  # Step #2
  ###############################################################################################
  
  ###################################################################################################
  #2. Step 2 indexes, i =1,2,..., I, each equivalent form of eta(i.e.,eta_i)  all.perm_index.
  # via a permutation matrix, P (i.e., Pi) 
  ###################################################################################################
  
  all.perm_index=as.matrix(lapply(1:length(s_n.factor), function(i) {as.matrix(cbind(matrix(s_n.factor[[i]], nrow=nrow(v_sign), ncol=n.factor, byrow=TRUE), v_sign))}))
  
  all.perm_index=cbind(1:n.perm, do.call(rbind, all.perm_index)) # Matrix to index each individual permutation. 
  
  colnames(all.perm_index)[1]="permutation index"
  names1=paste("order", 1:n.factor, sep="")
  colnames(all.perm_index)[2:(n.factor+1)]=names1
  
  #  write.table(all.perm_index, sep=",", file=paste("all.perm_index_", filename,".csv",sep=""), row.names=FALSE)  
  
  P=lapply(1:nrow(all.perm_index), function(i) {as(as.integer(c(all.perm_index[i,2:(n.factor+1)])), "pMatrix")%*%diag(all.perm_index[i,(n.factor+2):(2*n.factor+1)])}) 
  
  pop_data_lambda1=pop_lambda[1:n.var, 1:n.factor]
  pop_perm=lapply(1:length(P), function(i) {pop_data_lambda1%*%solve(P[[i]])}) # Create population 
  
  ###################################################################################################  
  # 3. Step 3 determines which eta_i each replication follows based on bias (first) and MSE (if ties were found based on bias) value. 
  # "correct_perm" shows which permutation index each replication should follow.
  ###############################################################################################
  sample_data_lambda1=list()
  sample_data_lambda1=lapply(1:length(sample_lambda), function(i) {sample_lambda[[i]][1:n.var, 1:n.factor]})
  
  bias_MSE_K_lambda=lapply(1:length(sample_data_lambda1), function(i) {cbind(matrix(seq(1:length(pop_perm))), do.call(rbind,lapply(1:length(pop_perm), function(j) {cbind(sum(sample_data_lambda1[[i]]-pop_perm[[j]])/(nrow(sample_data_lambda1[[i]])*ncol(sample_data_lambda1[[j]])), (tr(t(sample_data_lambda1[[i]]-pop_perm[[j]])%*%(sample_data_lambda1[[i]]-pop_perm[[j]])))/(nrow(sample_data_lambda1[[i]])*ncol(sample_data_lambda1[[i]])), sum(colSums(sample_data_lambda1[[i]]*pop_perm[[j]])/sqrt(colSums(sample_data_lambda1[[i]]^2)*colSums(pop_perm[[j]]^2)))/n.factor)})))})
  
  
  correct_perm=do.call(rbind,lapply(1:length(sample_data_lambda1), function(i) {ifelse(length(which(bias_MSE_K_lambda[[i]][,3]==min(bias_MSE_K_lambda[[i]][,3])))==1, which(bias_MSE_K_lambda[[i]][,3]==min(bias_MSE_K_lambda[[i]][,3])), which(abs(bias_MSE_K_lambda[[i]][,4])==min(abs(bias_MSE_K_lambda[[i]][,4]))))}))
  
  correct_perm=do.call(rbind,lapply(1:length(sample_data_lambda1), function(i) {cbind(i, ifelse(length(which(bias_MSE_K_lambda[[i]][,3]==min(bias_MSE_K_lambda[[i]][,3])))==1, which(bias_MSE_K_lambda[[i]][,3]==min(bias_MSE_K_lambda[[i]][,3])), which(abs(bias_MSE_K_lambda[[i]][,4])==min(abs(bias_MSE_K_lambda[[i]][,4])))))}))

  replication=matrix(correct_perm[,1])
  permutation=matrix(correct_perm[,2])
  f.rep.perm=table(replication, permutation)
  s.rep.perm=margin.table(table(replication, permutation), 2)

  correct_perm_list=do.call(rbind, lapply(1:length(permutation), function(i) {cbind(i, matrix(all.perm_index[permutation[i,], ], nrow=1), matrix(bias_MSE_K_lambda[[i]][permutation[i,],2:4], nrow=1))}))

  colnames(correct_perm_list)=c("rep#", "correct.perm", paste("order", 1:n.factor, sep=""), paste("sign", 1:n.factor, sep=""), "Bias", "MSE", "K")
  
  write.table(correct_perm_list, sep=",", file=paste("correct_perm_list.csv",sep=""), row.names=FALSE)
  
  ###############################################################################################
  # 4. Step 4: reorders and/or reflects the estimated matrix based on the corresponding a P matrix. Then, below will give correctly-computed bias & MSE based on the corresponding a P matrix. 
  ###############################################################################################
  
  pp=list()
  
  
  for (kk in 1:length(sample_lambda)){
    pp[[kk]]=diag(n.eta)
    pp[[kk]][1:n.factor, 1:n.factor]=P[[permutation[[kk]]]]
  }
  
  P=pp
  
  P=do.call(rbind,lapply(1:length(P), function(i) {P[[i]]}))
  
  output_list=list("n.perm"=n.perm, "permutation" =all.perm_index,"replication.permutation" = f.rep.perm, "summary.permutation"=s.rep.perm, "correct.permutation"=correct_perm_list)
  write.table(P, sep=",", file="P.csv", row.names=FALSE)
  return(output_list)
} 

