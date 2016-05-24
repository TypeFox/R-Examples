
genmvnorm<-function(cor,k,n,seed=F){
  
  
  if(is.matrix(cor)==F){
    x<-length(cor)
    if(x != (k*(k-1)/2) ){stop("STOP: wrong correlation table")}
    
    cr.cor<-matrix(NA,k,k)
    diag(cr.cor)<-1
    cr.cor[lower.tri(cr.cor)]<-cor
    cr.cor[upper.tri(cr.cor)]<-t(cr.cor)[upper.tri(cr.cor)]
    
    e<-eigen(cr.cor)
    L<-e$values #placing the eigenvalues in V
    Vm<-matrix(0,nrow=k,ncol=k) #creating a k x k matrix.
    diag(Vm)<-L #putting the eigenvalues on the diagonals
    Vm #check-- matrix with eigenvalues
    e$vectors #these are the eigenvectors
    l<-e$vectors %*% sqrt(Vm) #these are the loadings
  }
  
  if(is.matrix(cor)==T){ #if a correlation matrix was used
    
    e<-eigen(cor)
    L<-e$values #placing the eigenvalues in V
    Vm<-matrix(0,nrow=k,ncol=k) #creating a k x k matrix.
    diag(Vm)<-L #putting the eigenvalues on the diagonals
    Vm #check-- matrix with eigenvalues
    e$vectors #these are the eigenvectors
    l<-e$vectors %*% sqrt(Vm) #these are the loadings
    
  }
  
  if(seed != F){set.seed(seed)} 
  f<-matrix(nrow=k,l)
  #dim(f)
  ma<-matrix(nrow=n,ncol=k)
  for (i in 1:k){
    ma[,i]<-rnorm(n)
  }
  tma<-t(ma)
  sol<-f %*% tma
  sol<-t(sol)
  
  data<-data.frame(sol)
  return(data)
  
}

