getPerms<-function(n){
	mat<-matrix(0,2^n,n)
	for(i in 1:n)
		mat[,i]<-rep(0:1,times=2^(i-1),each=2^(n-i))
	mat
}

