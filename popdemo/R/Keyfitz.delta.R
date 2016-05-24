Keyfitz.delta <-
function(A,vector){
if(any(length(dim(A))!=2,dim(A)[1]!=dim(A)[2])) stop("A must be a square matrix")
if(!is.matrix_irreducible(A)) stop("Matrix A is reducible")
if(!is.matrix_primitive(A)) warning("Matrix A is imprimitive")
reigs<-eigen(A)
lmax<-which.max(Re(reigs$values))
w<-as.matrix(reigs$vectors[,lmax])
if(max(Im(w))>0) stop("Dominant right eigenvector contains nonzero imaginary components")
w<-abs(Re(w))
w<-w/sum(w)
vector<-as.matrix(vector/sum(vector))
delta<-0.5*norm(vector-w)
return(delta)
}

