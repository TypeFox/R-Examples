supermat <-
function(X){

size=dim(X)
n=size[1]
m=size[2]
p=size[3]
sup.mat.A=matrix(,nrow=n,ncol=m*p)
k=1
while (k<=p){
	 sup.mat.A[,((k-1)*m+1):(k*m)]=X[,,k]
	 k=k+1
	}
sup.mat.B=matrix(,nrow=m,ncol=p*n)
i=1
for (i in 1:n){
	sup.mat.B[,((i-1)*p+1):(i*p)]=X[i,,]
	}
sup.mat.C=matrix(,nrow=p,ncol=n*m)
j=1
while (j<=m){
	sup.mat.C[,((j-1)*n+1):(j*n)]=t(X[,j,])
	j=j+1
	}
out=list()
out$Xa=sup.mat.A
out$Xb=sup.mat.B
out$Xc=sup.mat.C
return(out)
}
