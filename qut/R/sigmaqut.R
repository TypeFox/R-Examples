sigmaqut <-
function(y,X,estimator='unbiased',intercept=TRUE,alpha.level='default',M=1000,qut.standardize=TRUE,penalty.factor=rep(1,p),offset=NULL,...){
#FUNCTION TO ESTIMATE SIGMA USING qut
	#Hidden function to obtain the cost for each sigma
	costfunction=function(sigma, fit1, fit2,  X1, X2, y1, y2, lambda1, lambda2, N, estimator,O1,O2){

		#Variable selection for each set
		betahat1lasso=coef(fit1, sigma*lambda1/N,offset=O1)[-1]
		betahat2lasso=coef(fit2, sigma*lambda2/N,offset=O2)[-1]
		w1=which(betahat1lasso!=0)
		w2=which(betahat2lasso!=0)
		
		out=NULL
		out$w1=NA
		out$w2=NA
		if( (length(w1)>N) || (length(w2)>N)) cf=Inf
		else{
			if(length(w1)>=1) fit1OLS=lm(y2~X2[,w1],offset=O1) else fit1OLS=lm(y2~1,offset=O1)
			if(length(w2)>=1) fit2OLS=lm(y1~X1[,w2],offset=O2) else fit2OLS=lm(y1~1,offset=O2)
			if(estimator=='unbiased'){ #Estimate sigma with unbiased estimator
				
				sigma2.1= sum(fit1OLS$res^2)/(length(s1)-length(w1)-1)
				sigma2.2= sum(fit2OLS$res^2)/(length(s1)-length(w2)-1)
			}
			else if(estimator=='mle'){ #Estimate sigma with maximum likelihood estimator
				sigma2.1=sum(fit1OLS$res^2)/(length(s1))
				sigma2.2=sum(fit2OLS$res^2)/(length(s1))
			}
	   
			#sigmahat is the mean of the 2 estimators
			sigmahat=sqrt(0.5*(sigma2.1+sigma2.2))
			cf=abs(sigmahat-sigma)

			if(sigma2.1==0 || sigma2.2==0 || is.na(sigma2.1) || is.na(sigma2.2)) cf=Inf
			out$w1=length(w1)
			out$w2=length(w2)
		}
	
		out$cf=cf
		return(out)
	}

	######
	
	p=ncol(X);P=p
	n=nrow(X);N=n

	no.penalty=which(penalty.factor==0)
	if(length(no.penalty)==0) no.penalty=NULL

	#Check for warnings
	if (is.null(p) | (p <= 1)) stop("X should be a matrix with 2 or more columns")
	if (n!=length(y)) stop("Number of observations in y not equal to number of rows in X")

	#initialize sigmas
    fit0=lm(y~1,offset=offset)
	sigmamax=sqrt(mean(fit0$res^2))
	ntry=30
	sigmas=seq(1.e-6, sigmamax, length=ntry)
	lNA=5
	sigmahats=rep(NA, lNA)
	
	#Estimate sigma lNA times to be robust
	for(i in 1:lNA){
		#Split data in 2
		s1=sample(n,n/2)
		X1=X[s1,]; y1=y[s1]
		X2=X[-s1,]; y2=y[-s1]
		
		if (is.null(offset)){O1=NULL ;O2=NULL}
		else{O1=offset[s1];O2=offset[-s1]}
		fit1=glmnet(X1,y1,family="gaussian",intercept=intercept,standardize=F,offset=O1,...)
		fit2=glmnet(X2,y2,family="gaussian",intercept=intercept,standardize=F,offset=O2,...)
		#lambda qut for each data 
		outqut=lambdaqut(y=y1,X=X1,alpha.level=alpha.level,M=M,qut.standardize=qut.standardize,family=gaussian,intercept=intercept,no.penalty=no.penalty,offset=O1)
		lambdaqut1=outqut$lambda
		X1=outqut$Xnew
		
		outqut=lambdaqut(y=y2,X=X2,alpha.level=alpha.level,M=M,qut.standardize=qut.standardize,family=gaussian,intercept=intercept,no.penalty=no.penalty,offset=O2)
		lambdaqut2=outqut$lambda
		X2=outqut$Xnew
		
		#Regularization path for each data
		fs=rep(NA, length(sigmas))
		jopt=NA;fsmin=Inf
		
		#maximize cost function for a grid of sigmas
		for(j in 1:ntry){
			temp=costfunction(sigmas[j], fit1=fit1, fit2=fit2, X1=X1, X2=X2, y1=y1, y2=y2, lambda1=lambdaqut1,lambda2=lambdaqut2, N=N/2,estimator=estimator,O1=O1,O2=O2)
			if(temp$cf<fsmin){
				fsmin=temp$cf
				jopt=j
			}
			fs[j]=temp$cf
		}
		if(is.na(jopt)){
			warning("No sigma was able to be estimated, check for data problems")
			sigmahat=NA
		}
	    else sigmahat=sigmas[jopt]
		sigmahats[i]=sigmahat
	}

	#robust estimator
	sigmahat=median(sigmahats)

	out=sigmahat
	return(out)
}
