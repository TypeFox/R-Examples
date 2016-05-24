orth <-
function(A){

svd=svd(A)
m=nrow(A)
n=ncol(A)
if (m>1){
	s=svd$d
} else{ 
	if (m==1){
		s=svd$d[1]
		} else{
		s==0
	}
}
eps=1.0e-016
tol=max(m,n)*max(s)*eps
r=sum(s>tol)
Q=svd$u[,1:r]
return(Q)
}
