sigmarcv <-
function(y,X,cv=FALSE,fit=NA,intercept=TRUE){

#Taken from:

#Variance estimation using refitted cross-validation in
#ultrahigh dimensional regression

#Jianqing Fan,
#Princeton University, USA
#Shaojun Guo
#Chinese Academy of Sciences, Beijing, PRC
#and Ning Hao
#University of Arizona, Tucson, USA

#Journal of Royal Statistics Society
#2011


#kindly proportioned by the paper authors
 

	n=nrow(X)
	p=ncol(X)
    
	# we assume the data have been resampled by random permutation. If not, resample the label here.

	k          <- floor(n/2)
    x1         <- X[1:k,]
    y1         <- y[1:k]
    x2         <- X[(k+1):n,]
    y2         <- y[(k+1):n]
	
	out=NULL
	 
	if(cv==TRUE){
		#CV
		flasso=fit
		sigma.est1=NA
		maxiter=0
		while(is.na(sigma.est1)){
			if(is.na(fit)) flasso     <- lars(X,y,type="lasso",eps=1e-10,use.Gram=F,intercept=intercept);
			object     <- cv.lars(X,y,type="lasso",use.Gram=F,plot.it=F,se=FALSE)
			beta.est   <- predict.lars(flasso,s=object$index[which.min(object$cv)],
							 type="coefficients",mode="fraction")$coef

			df         <- length(which(abs(beta.est)>0))
			if(intercept) beta0=mean(y)-mean(X%*%beta.est) else beta0=0

			sigma.est1 <- sum((y-X%*%beta.est-beta0)^2)/(n-df)
			out$sigmahat=sqrt(sigma.est1)
			out$ind=which(beta.est!=0)
			maxiter=maxiter+1
			if(maxiter==10) stop("The maximum number of iterations where reached and all cross validation models were full rank. No estimation of sigma was able to be done by sigma='cv'.  Try sigma='qut'.")
		}
	}
	else if(cv==FALSE){
		#RCV
		maxiter=0
		sigma.est3=NA
		sigma.est4=NA
		while(is.na(sigma.est3)|is.na(sigma.est4)){
			flasso     <- lars(x1,y1,type="lasso",eps=1e-10,use.Gram=F,intercept=intercept);
			object     <- cv.lars(x1,y1,type="lasso",use.Gram=F,plot.it=F,se=FALSE)
			beta.est1  <- predict.lars(flasso,s=object$index[which.min(object$cv)],
						 type="coefficients",mode="fraction")$coef
			ind1       <- which(abs(beta.est1)>0)
			if ( length(ind1) == 0 ) {
				sigma.est3 <- var(y2)
			} else {
				if(intercept) object     <- lm(y2 ~ x2[,ind1] + 1) else object     <- lm(y2 ~ x2[,ind1] - 1)
				sigma.est3 <- sum((object$resid)^2)/(n-k-length(ind1))
			}

			flasso     <- lars(x2,y2,type="lasso",eps=1e-10,use.Gram=F,intercept=intercept);
			object     <- cv.lars(x2,y2,type="lasso",use.Gram=F,plot.it=F,se=FALSE)
			beta.est2  <- predict.lars(flasso,s=object$index[which.min(object$cv)],
							type="coefficients",mode="fraction")$coef
			ind2       <- which(abs(beta.est2)>0)
			if ( length(ind2) ==0 ) {
				sigma.est4 <- var(y1)
			} else {
				if(intercept) object     <- lm(y1 ~ x1[,ind2] + 1) else object     <- lm(y1 ~ x1[,ind2] - 1) 
				sigma.est4 <- sum((object$resid)^2)/(k-length(ind2))
			}
			maxiter=maxiter+1
			if(maxiter==10) stop("The maximum number of iterations where reached and all RCV models were full rank. No estimation of sigma was able to be done by sigma='rcv'.  Try sigma='qut'.")
		}
		out$sigmahat=sqrt((sigma.est3+sigma.est4)/2)
		out$ind=which(abs(beta.est1*beta.est2)>0)
	}
	
	return(out)
 }