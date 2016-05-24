showPCA <- function(result,n_joint=0,n_indiv=rep(0,length(result$data)),Colors='black'){
  
  old.par <- par(no.readonly = TRUE) # all par settings which could be changed
  on.exit(par(old.par))
  l <- length(result$data)
  
  nPCs = n_joint+sum(n_indiv)
  PCs = matrix(nrow=nPCs,ncol = dim(result$data[[1]])[2])
  PC_names = rep('',nPCs)
  if(n_joint>0){
    SVD = svd(do.call(rbind,result$joint),nu=n_joint,nv=n_joint)
    PCs[1:n_joint,] = diag(SVD$d)[1:n_joint,1:n_joint]%*%t(SVD$v[,1:n_joint])
    PC_names[1:n_joint] = paste("Joint ",1:n_joint) 
  }
  
  for(i in 1:l){
  if(n_indiv[i]>0){
    SVD = svd(result$individual[[i]],nu=n_indiv[i],nv=n_indiv[i])
    indices = (n_joint+sum(n_indiv[0:(i-1)])+1):(n_joint+sum(n_indiv[0:i]))
    PCs[indices,] = diag(SVD$d)[1:n_indiv[i],1:n_indiv[i]]%*%t(SVD$v[,1:n_indiv[i]])
    PC_names[indices] = paste(names(result$data)[i]," Indiv ",1:n_indiv[i]) 
  }}
  nplots = (nPCs-1)^2
  par(mar = c(4,4,2,2))
  layout(matrix(c(1:nplots),nrow = nPCs-1,ncol = nPCs-1))
  for(i in 2:nPCs){
    for(j in 1:(nPCs-1)){
    if(j>=i) plot.new()
    else plot(PCs[i,],PCs[j,],xlab = PC_names[i],ylab = PC_names[j],col=Colors)
  }}
}

