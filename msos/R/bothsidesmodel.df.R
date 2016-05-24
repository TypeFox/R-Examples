bothsidesmodel.df <-
function(xx,n,pattern) {
	l <- ncol(pattern)
	tt <- pattern==1
	pj <- apply(pattern,2,sum)
	
	df <- matrix(0,l,l)
	diag(df) <- n-pj
	
	a <- vector("list",l)
	
	a[[1]] <- solve(xx[tt[,1],tt[,1]])
	
	for(i in 2:l) {
		if(pj[i]>0) a[[i]] <- solve(xx[tt[,i],tt[,i]])
		for(j in 1:(i-1)) {
			if(pj[i]==0|pj[j]==0) {pij<-0}
			else {
				b <- xx[tt[,i],tt[,j]]
				if(is.vector(b)) b <- matrix(b,nrow=pj[i])
				pij <- tr(a[[i]]%*%b%*%a[[j]]%*%t(b))
			}
			df[j,i] <- df[i,j] <- n-pj[i]-pj[j]+pij
		}
	}
	df
}
