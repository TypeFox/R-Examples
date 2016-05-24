`FromSymmetricStorageUpper` <-
function(x){
 n<-floor((-1+sqrt(1+8*length(x)))/2)
 z<-matrix(numeric(n^2), nrow=n)
 i<-as.vector(lower.tri(z,diag=TRUE))
 z[i]<-x
 ztranspose<-t(z)
 diag(ztranspose)<-0
 z+ztranspose
}

