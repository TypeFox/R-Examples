multCh7<-function(our.matrix){
  k<-ncol(our.matrix)
  n<-nrow(our.matrix)
  n.perm<-factorial(k)^n
  outp<-array(dim=c(n,k,n.perm))
  sorted.rows<-t(apply(our.matrix,1,sort))
  possible.row.arr<-array(dim=c(factorial(k),k,n))
  for(i in 1:n){
    possible.row.arr[,,i]<-matrix(unlist(permn(sorted.rows[i,])),byrow=T,ncol=k)
  }  
  
  get.mat<-function(index){
    possible.mat<-NULL
    for(j in 1:n){
      possible.mat<-rbind(possible.mat,possible.row.arr[index[j],,j])
    }
    return(possible.mat)	
  }
  
  
  our.list<-list(1:factorial(k))
  for(i in 2:n){
    our.list[[i]]<-1:factorial(k)
  }
  index.grid<-as.matrix(expand.grid(our.list),ncol=n)
  
  for(i in 1:n.perm){
    outp[,,i]<-get.mat(index.grid[i,])
  }
  outp
}