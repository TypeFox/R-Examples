simul.mp<-function(n,m,mu=0,sigmap=1,p=2){
ris<-list()
a<-rnormp(n*m,mu=mu,sigmap=sigmap,p=p)
 a<-matrix(a,n,m)
 xc<-split(a,col(a))
  mat<-matrix(nrow=m,ncol=6)
   for (i in 1:m){
    A<-paramp(xc[[i]])
    for (j in 1:6) mat[i,j]<-A[[j]]
}
ma<-matrix(nrow=2,ncol=5,dimnames=list(c("Mean","Variance"),c("Mean","Mp","Sd","Sp","p")))
  for(j in 1:5){ma[1,j]<-mean(mat[,j]);ma[2,j]<-(length(mat[,j])-1)*var(mat[,j])/length(mat[,j])}
ris$dat<-mat
ris$table<-ma
ris$iter<-sum(mat[,6])
class(ris)<-"simul.mp"
ris
}


print.simul.mp<-function(x,...){
print(x$table)
cat("\nNumber of samples with a difficult convergence:",x$iter,"\n")
invisible(x)
}

