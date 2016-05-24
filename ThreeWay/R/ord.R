ord <-
function(X){

if (is.vector(X)==TRUE){
	A=sort(X,method="sh",index.return=TRUE)$x
	a=sort(X,method="sh",index.return=TRUE)$ix
}
if (is.vector(X)==FALSE){
	A=matrix(,nrow(X),ncol(X))
	a=matrix(,nrow(X),ncol(X))
	for (j in 1:ncol(X)){ 
		A[,j]=sort(X[,j],method="sh",index.return=TRUE)$x
		a[,j]=sort(X[,j],method="sh",index.return=TRUE)$ix
	}
}
out=list()
out$A=A
out$a=a
return(out)
}
