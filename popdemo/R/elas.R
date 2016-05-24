elas <-
function(A,eval="max"){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(eval=="max"){
    val<-which.max(abs(Re(eigen(A)$values)))
}
else{
    val<-eval
}
reigs<-eigen(A)
leigs<-eigen(t(A))
lambda<-reigs$values[val]
w<-as.matrix(reigs$vectors[,val])
v<-as.matrix(leigs$vectors[,val])
S<-(v%*%t(w))/as.vector(t(v)%*%w)
E<-(1/lambda)*S*A
if(max(Im(E))>0) return(E) else(return(Re(E)))
}

