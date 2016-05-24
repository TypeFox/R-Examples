proportion.function<-function(k, beta.clu.ward, beta.clu.complete, dist, return.list=F){
  k.w=k[1]
  k.c=k[2]
  tree.ward<-cutree(beta.clu.ward,k=k.w)
  tree.complete<-cutree(beta.clu.complete,k=k.c)
  N<-sum(tree.ward*0+1)  
  A<-table(tree.ward, tree.complete)
  A.star<-diag(0,k.w)
  if(return.list==T) beta.list<-list()
  i=1
  for(i in 1:k.w){
    A.star[i,i]<-max(A)
    A.max<-which(A == max(A), arr.ind = TRUE)[1,]
    r.max<-A.max[1]
    c.max<-A.max[2]
    if(return.list==T) beta.list[[i]]<-(1:N)[(tree.ward==r.max)&(tree.complete==c.max)]
    A[r.max,]<-0
    A[,c.max]<-0
  } 
  while ((A.star[nrow(A.star), ncol(A.star)]==0)) {
    A.star=A.star[-(nrow(A.star)), -(ncol(A.star))]
    if(return.list==T) 
      beta.list = beta.list[-length(beta.list)]
    if(is.null(dim(A.star)))
      break
  }
  if(return.list==T) return(list("beta.list"=beta.list, "A.star"=A.star))
  sum(A.star)
}
