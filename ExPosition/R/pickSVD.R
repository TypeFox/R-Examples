pickSVD <-
function(datain,is.mds=FALSE,decomp.approach='svd',k=0){

	dataDims <- dim(datain)
	I <- dataDims[1]
	J <- dataDims[2]
	m <- min(I,J)

	#check k
	if(k < 1 || k > m){
		k <- m
	}
	
	
	flip <- FALSE	
	if (I < J){
		datain <- t(datain)
		flip <- TRUE
	}	
	
	#check decomp.approach
	if(is.null(decomp.approach)){
		decomp.approach <- 'svd'
	}
	#if there are more than 1 million elements or you want to go fast, do the eigen decomp.approach!
	num.el <- (I * J)
	if(num.el > 1000000){
		decomp.approach <- 'eigen'
	}
	
	##for now, only eigen & svd.
	if(tolower(decomp.approach)=='eigen'){
		eigOut <- eigen(t(datain) %*% datain)	
		Q <- eigOut$vectors
		d <- sqrt(eigOut$values)
		P <- datain %*% Q %*% diag(d^-1)

	}else{ ##the default method.
		svd.out <- svd(datain,nu=k,nv=k)
		P <- svd.out$u
		Q <- svd.out$v
		d <- svd.out$d
	}
	##but we hope to add faster methods soon, e.g., RcppArmadillo
	
	
	if(flip){
		temp<-Q
		Q<-P
		P<-temp
		rownames(Q) <- rownames(datain)
		rownames(P) <- colnames(datain)		
	}else{
		rownames(P) <- rownames(datain)
		rownames(Q) <- colnames(datain)			
	}
	
	#this guarantees I take the rank as determined by the sings/eigs/"tau"
	P <- P[,1:min(length(d),ncol(P))]
	Q <- Q[,1:min(length(d),ncol(Q))]
		
	#find precision limit, fix what comes back.
	precisionLimit <- 2*.Machine$double.eps	
	
	if(is.mds){
		indToKeep <- which(d > precisionLimit)	
		tau <- d[indToKeep]/sum(d[indToKeep])	##value could be small due to error.
	}else{
		indToKeep <- which(d^2 > precisionLimit)		
		tau <- d[indToKeep]^2/sum(d[indToKeep]^2)	##value could be small due to error.
	}

	indToKeep <- indToKeep[1:min(c(length(indToKeep),k))]	
	return(list(u=P[,indToKeep],v=Q[,indToKeep],d=d[indToKeep],tau=tau*100))
}
