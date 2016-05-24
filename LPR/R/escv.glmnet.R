escv.glmnet<-function(x, y, lambda, nfolds=10, foldid, tau=0, mLS=FALSE, nu=0, 
			standardize=TRUE, intercept=TRUE,...){
# Inputs:
#  x,y: data set
#  lambda: tuning parameter of Lasso problem (argmin ||Y-X\beta||^2+lambda*||\beta||_1.
#  nfolds: number of folds.
#  foldid: the id of different folds
#  mLS: if mLS=TRUE, I will use Lasso+mLS estimator. The default value is FALSE
#  tau: the tuning parameter in modified Least Squares (mLS)
#  nu: the tuning parameter in Ridge
#  standardize,intercept,...:  parameters for glmnet

# Outputs:
#  cv: cv values
#  cv.error: standard error of cv
#  es: es values
#  es.error: standard error of esnfolds
#  lambda: tuning parameter
#  lambda.cv: the lambda selected by cv
#  lambda.cv1se: the lambda selected by cv1se
#  lambda.escv: the lambda selected by escv

	n<-dim(x)[1]
	p<-dim(x)[2]
	if(missing(lambda)){
		lambda_max<-max(abs(t(x)%*%y)/n)
		llam<-seq(log(lambda_max),log(lambda_max*0.0001),length=100)
		lambda<-exp(llam)
	}	
	
	if(missing(foldid))
		foldid <- sample(rep(seq(nfolds), length = n))
	else nfolds <- max(foldid)

	out<-list()
		for(k in 1:nfolds){
			test<-foldid==k
			train<-foldid!=k
			obj<-glmnet(x[train,,drop=FALSE], y[train], lambda=lambda, standardize=standardize, intercept=intercept,...)
			fitmat<-predict(obj,newx=x)
			predtest<-predict(obj,newx=x[test,,drop=FALSE])				
#			if(mLS | Ridge){
#				select.all<-obj$beta!=0
#			} else {
#				select.all<-NULL
#			}	
			residmat<-apply((y[test] - predtest)^2, 2, mean)
			out[[k]]<-list(residmat=residmat, fitmat=fitmat)			
#			out[[k]]<-list(residmat=residmat, fitmat=fitmat, select.all=select.all)
		}			

	
	residmat<-matrix(0,length(lambda),nfolds)
	residmates<-matrix(0,length(lambda),nfolds)
	fitmat<-array(0,dim=c(n,length(lambda),nfolds))
		
	for(k in 1:nfolds){
		residmat[,k]<-out[[k]]$residmat
		fitmat[,,k]<-out[[k]]$fitmat
		out[[k]]$residmat<-0
		out[[k]]$fitmat<-0
	}

	meanfit<-apply(fitmat,c(1,2),mean)
	meanfit2<-apply(meanfit^2,2,sum)
	for(k in 1:nfolds){
		residmates[,k]<-apply((fitmat[,,k]-meanfit)^2,2,sum)/meanfit2
	}
	residmates[is.na(residmates)]<-Inf

	cv<-apply(residmat,1,mean)
	cv.error<-sqrt(apply(residmat,1,var)/nfolds)	
	es<-apply(residmates,1,mean)
	es.error<-sqrt(apply(residmates,1,var)/nfolds)
	
	indcv<-which.min(cv)
	lambda.cv<-lambda[indcv]
	
	cv1se<-cv
	cv1se[cv<=(cv[indcv]+cv.error[indcv])]<-cv[indcv]+cv.error[indcv]
	indcv1se<-which.min(cv1se)
	lambda.cv1se<-lambda[indcv1se]	
	
	indescv<-which.min(es[1:indcv])
	lambda.escv<-lambda[indescv]
	
	# if mLS or Ridge, refitting the model by mLS or Ridge
	if(mLS){
		out<-list()
			for(k in 1:nfolds){
				test<-foldid==k
				train<-foldid!=k
				obj<-glmnet(x[train,,drop=FALSE], y[train], lambda=lambda[1:indcv], 
					standardize=standardize, intercept=intercept,...)
				fitmat<-predict(obj,newx=x)
				predtest<-predict(obj,newx=x[test,,drop=FALSE])
				if(mLS){  # if mLS=TRUE, I will use Lasso to select predictors and then use mLS to refit
					selectset0<-rep(0,p)
					for(i in 1:indcv){
						selectset<-abs(obj$beta[,i])>0
						if(sum(selectset)>0){
							if(sum(abs(selectset-selectset0))>0){
								mls.obj<-mls(x[train,selectset,drop=FALSE], y[train], tau=tau,
									standardize=standardize, intercept=intercept)
								fitmat[,i]<-drop(scale(x[,selectset,drop=FALSE], mls.obj$meanx, FALSE) %*% t(mls.obj$beta)) + mls.obj$mu
								predtest[,i]<-drop(scale(x[test,selectset,drop=FALSE], mls.obj$meanx, FALSE) %*% t(mls.obj$beta)) + mls.obj$mu
							} else{
								fitmat[,i]<-fitmat[,i-1]
								predtest[,i]<-predtest[,i-1]
							}
							selectset0<-selectset
						}		
					}
				}
				residmat<-apply((y[test] - predtest)^2, 2, mean)
				out[[k]]<-list(residmat=residmat,fitmat=fitmat)
			}			
		
		
		residmat<-matrix(0,indcv,nfolds)
		residmates<-matrix(0,indcv,nfolds)
		fitmat<-array(0,dim=c(n,indcv,nfolds))
		
		for(k in 1:nfolds){
			residmat[,k]<-out[[k]]$residmat
			fitmat[,,k]<-out[[k]]$fitmat
			out[[k]]$residmat<-0
			out[[k]]$fitmat<-0
		}

		meanfit<-apply(fitmat,c(1,2),mean)
		meanfit2<-apply(meanfit^2,2,sum)
		for(k in 1:nfolds){
			residmates[,k]<-apply((fitmat[,,k]-meanfit)^2,2,sum)/meanfit2
		}
		residmates[is.na(residmates)]<-Inf

		cv<-apply(residmat,1,mean)
		cv.error<-sqrt(apply(residmat,1,var)/nfolds)	
		es<-apply(residmates,1,mean)
		es.error<-sqrt(apply(residmates,1,var)/nfolds)
	
		indcv<-which.min(cv)
		lambda.cv<-lambda[indcv]
	
		cv1se<-cv
		cv1se[cv<=(cv[indcv]+cv.error[indcv])]<-cv[indcv]+cv.error[indcv]
		indcv1se<-which.min(cv1se)
		lambda.cv1se<-lambda[indcv1se]	
	
		indescv<-which.min(es[1:indcv])
		lambda.escv<-lambda[indescv]
						
	}
	
	object<-list(lambda=lambda, cv=cv, cv.error=cv.error, es=es, es.error=es.error, lambda.cv=lambda.cv,
		lambda.cv1se=lambda.cv1se, lambda.escv=lambda.escv)
	invisible(object)	
}