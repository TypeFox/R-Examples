np.t=function(a,b,c,d,n,alpha){
		if (a+b+c+d!=n){return(paste("Caution: a+b+c+d is not equal to n."))}	
		if (a+b+c+d==n){
		data.dich=matrix(c(rep(1,(a+b)),rep(0,(c+d)),rep(1,a),rep(0,b),rep(1,c),rep(0,d)),ncol=2, nrow=n)
		ntotal=2*n
		X0=matrix(ncol=2,nrow=n,0.5)
		U1sens=data.dich[,1]
		U2sens=data.dich[,2]
		data.dich=rbind(data.dich,X0)
		ranks=matrix(ncol=2,nrow=ntotal,NA)
		ranks[,1]=rank(data.dich[,1],ties.method="average")
		ranks[,2]=rank(data.dich[,2],ties.method="average")
		ranks1=ranks[1:n,]
		ranks0=ranks[(n+1):ntotal,]
		value=1/ntotal*(apply(ranks1,2,mean)-apply(ranks0,2,mean))+1/2 	
		seU1=value[1]
		seU2=value[2]
		est=seU1-seU2
		ranks11=matrix(ncol=2,nrow=n,NA)
		ranks11[,1]=rank(data.dich[1:n,1],ties.method="average")
		ranks11[,2]=rank(data.dich[1:n,2],ties.method="average")
		ranks00=matrix(ncol=2,nrow=n,(n+1)/2)
		z0 = ranks0-ranks00
		z1 = ranks1-ranks11
		z0quer=apply(ranks0,2,mean)-apply(ranks00,2,mean)
		z1quer=apply(ranks1,2,mean)-apply(ranks11,2,mean)
		zsum0 <- matrix(0,2,2)
		zsum1 <- matrix(0,2,2)
		for (j in 1:n){
			vi0=(z0[j,]-z0quer)%*%t(z0[j,]-z0quer)
			zsum0 <- zsum0+vi0
		}

		for (j in 1:n){
			vi1=(z1[j,]-z1quer)%*%t(z1[j,]-z1quer)
			zsum1 <- zsum1+vi1
		}
		v0=ntotal/((ntotal-n)^2*n*(n-1))*zsum0
		v1=ntotal/((ntotal-n)^2*n*(n-1))*zsum1
		v=v0+v1
		contr=c(1,-1)
		
		df=max(1,sum(diag(v))^2/(v[1,1]^2/(n-1)+v[2,2]^2/(n-1)))
		if (df=="NaN"){c(0,0)}
		if (df!="NaN"){
		np_t=value%*%contr+c(-1,1)*(sqrt(contr%*%v%*%as.matrix(contr,ncol=1,nrow=2))*qt(1-alpha/2,df=df))/sqrt(ntotal)}
		attr(np_t, "conf.level") <- 1-alpha
		rval <- list(conf.int = np_t, estimate = est)
		class(rval) <- "htest"
		return(rval)
}}