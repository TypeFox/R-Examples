nrm2 <-
function(A){

n=nrow(A)
m=ncol(A)
d=SUM(A)$col
N=matrix(0,n,m)
if (sum(d)>0){
	for (k in 1:m){
		if (d[k]>1e-30){
			N[,k]=A[,k]/(matrix(1,n)*d[k]^.5)
		}
	}
}
return(N)
}
