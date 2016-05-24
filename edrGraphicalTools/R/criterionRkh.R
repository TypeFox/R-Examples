criterionRkh <-
function(Y,X,H,K,indices,B=50,method){

	cat("Be patient. The program is running ... \n")
	cl <- match.call()
	args <- as.list(match.call(criterionRkh))[-1]

	if((any(is.na(X))==T)|(any(is.na(Y))==T))
		stop("The data should not contained missing value")
	if(length(Y)!=dim(X)[1]) 
		stop("Y and X should have the same number of row")
	n<-length(Y)
	p<-ncol(X)
	if(missing(H)) 
		H<-round(seq(2,n/4,l=8))
	if(missing(K)) 
		K <- 1:(min(p,25))
	if(missing(method))
		stop("the method should be specified")
	if(!(method %in% c("SIR-I","SIR-II","SAVE")))
		stop("the method should be specified by SIR-I or SIR-II or SAVE")
	nH<-length(H)
	sigma<-((n-1)/n)*var(X)

	#generating bootstrap sample
	if (missing(indices)){
		indices<-sample(n,B*n,replace=T)
	} else { 
		if(length(indices)!=(B*n)){ 
			indices<-sample(n,B*n,replace=T)
			cat("The bootstrap sample is generating", "\n")
	} }

	projector <-function(beta1,sigm){
		beta1%*%(solve(t(beta1)%*%sigm%*%beta1))%*%t(beta1)%*%sigm
	}

	#Rest<-matrix(0,ncol=p,nrow=B)
	#resRkhboot<-matrix(0,ncol=p,nrow=nH)
	nK <- length(K)
	Rest<-matrix(0,ncol=nK,nrow=B)
	resRkhboot<-matrix(0,ncol=nK,nrow=nH)

	simboot<-list()
	i<-0

	###loop in H
	for (h in H){
		i<-i+1
		result<-edr(Y,X,H=h,K=1,method=method)

		###loop in sample Bootstrap
		for (iboot in 1:B) {
			Yboot<-Y[indices[c((n*(iboot-1)+1):(n*iboot))]]
			Xboot<-X[indices[c((n*(iboot-1)+1):(n*iboot))],]
			sigmaboot<-((n-1)/n)*var(Xboot)
			result_boot<-edr(Yboot,Xboot,H=h,K=1,method=method)

			#estimation of Rkh by bootstrap
			for (k in K){
				pkest<-projector(result$matEDR[,1:k],sigma)
				pkestboot<-projector(result_boot$matEDR[,1:k],sigmaboot)
				Rest[iboot,k]<-Re(sum(diag(pkest%*%pkestboot))/k)
		}	}

		simboot<-rbind(simboot,list(Rest))
		resRkhboot[i,]<-matrix(apply(Rest,2,mean),nrow=1)
	}

	res <- list(Rkhbootmean=Re(resRkhboot), Rkhboot=simboot, method=method, n=n, H=H,
							K=K, indices=indices)
	class(res) <-"criterionRkh"  
	res
}

#Debuggage de "proejctor"
	#projector <-function(beta1,sigm){
		#result <- try(beta1%*%(solve(t(beta1)%*%sigm%*%beta1))%*%t(beta1)%*%sigm)
		#if (inherits(result, "try-error")) {
		#	print(paste("Dimension of 'beta':", paste(dim(beta1), collapse=", ")))
		#	print(paste("Value of 'h' studied:",h))
		#	print("Matrix to invert:")
		#	print(t(beta1)%*%sigm%*%beta1)
		#	print("The covariance matrix has the following diagonal elements.")
		#	print(svd(sigm)$d)
		#	print(sliceMat(Yboot, Xboot, h))
		#	test1 <<- Yboot
		#	test2 <<- Xboot
		#	stop("Error when trying to invert a matrix.")
		#} else {
		#	result
		#}
	#}


