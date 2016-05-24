multCh7SM<-function(our.matrix){
  k.vals<-rowSums(!is.na(our.matrix))
  k<-max(k.vals)
  n<-nrow(our.matrix)
  n.perm<-prod(factorial(k.vals))
  outp<-array(dim=c(n,k,n.perm))
  sorted.rows<-t(apply(our.matrix,1,rank,na.last="keep"))
  
  possible.row.arr<-list()
  
  for(i in 1:n){
    possible.row.arr[[i]]<-matrix(unlist(permn(sorted.rows[i,!is.na(sorted.rows[i,])])),byrow=T,ncol=k.vals[i])
  }  
  
  get.mat<-function(index){
    possible.mat<-matrix(nrow=n,ncol=max(k.vals))
    for(j in 1:n){
      possible.mat[j,!is.na(our.matrix[j,])]<-possible.row.arr[[j]][index[j],]
    }
    return(possible.mat)	
  }
  
  
  our.list<-list(1:factorial(k.vals[1]))
  for(i in 2:n){
    our.list[[i]]<-1:factorial(k.vals[i])
  }
  index.grid<-as.matrix(expand.grid(our.list),ncol=n)
  
  for(i in 1:n.perm){
    outp[,,i]<-get.mat(index.grid[i,])
  }
  outp
}