Cohen.cumulative <-
function(A,vector){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
order<-dim(A)[1]
if(!is.matrix_irreducible(A)) stop("Matrix A is reducible")
if(!is.matrix_primitive(A)) warning("Matrix A is imprimitive")
M<-A
reigs<-eigen(M)
leigs<-eigen(t(M))
lmax<-which.max(Re(reigs$values))
lambda<-reigs$values[lmax]
A<-A/lambda
w<-as.matrix(reigs$vectors[,lmax])
v<-as.matrix(leigs$vectors[,lmax])
if(max(Im(w))>0|max(Im(v))>0) stop("Dominant eigenvectors contain nonzero imaginary components")
w<-abs(Re(w))
v<-abs(Re(v))
w<-w/sum(w)
v<-v/as.vector(t(v)%*%w)
vector<-vector/sum(vector)
I<-diag(order)
wv<-w%*%t(v)
D1v<-(solve(I+wv-A)-wv)%*%vector
D1<-sum(abs(D1v))
return(D1)
}

