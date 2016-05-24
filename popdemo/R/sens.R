sens <-
function(A,eval="max",all=FALSE){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
reigs<-eigen(A)
leigs<-eigen(t(A))
if(eval=="max"){
    val<-which.max(Re(reigs$values))
}
else{
    val<-eval
}
w<-as.matrix(reigs$vectors[,val])
v<-as.matrix(leigs$vectors[,val])
S<-(v%*%t(w))/as.vector(t(v)%*%w)
if(!all) S[A==0]<-0
if(max(Im(S))>0) return(S) else(return(Re(S)))
}

