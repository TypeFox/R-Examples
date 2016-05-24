qut <-
function(y,X,fit,family=gaussian,alpha.level='default',M=1000,qut.standardize=TRUE,intercept=TRUE,offset=NULL,bootstrap=TRUE,sigma=ifelse(n>2*p,'ols','qut'),beta0='iter',estimator='unbiased',type=c('glmnet','lars'),lambda.seq=0,penalty.factor=rep(1,p),lambda.min.ratio=ifelse(n<p,0.01,0.0001),nlambda=100,lambda=NULL,...){
#FUNCTION TO FIT GLM WITH THE QUANTILE UNIVERSAL THRESHOLD

	#Hidden function to get glm fit with non-zero coefficients found by LASSO-GLM
	betaGLMestimates=function(){
		if(type=='glmnet') betatemp=beta[-1] else betatemp=beta
		betaglm=betatemp*0
		converged=TRUE
		
		if(is.null(offset)) O=rep(0,n)
		if(sum(betatemp!=0)>0){
			if(intercept) outglm=glm(y~X0[,betatemp!=0],family=family,offset=offset)
			else outglm=glm(y~X0[,betatemp!=0]-1,family=family,offset=offset)
			betaglm0=outglm$coef
			converged=outglm$converged
		}
		else if(intercept){
			outglm=glm(y~1,family=family,offset=offset)
			betaglm0=outglm$coef
			converged=outglm$converged
		}
		
		if(intercept){
			betaglm[betatemp!=0]=betaglm0[-1]
			betaglm=c(betaglm0[1],betaglm)
		}
		else{
			betaglm[betatemp!=0]=betaglm0
			betaglm=c(0,betaglm)
		}
		
		names(betaglm)[1]='(Intercept)'	
		names(betaglm)[-1]=1:p
		attr(betaglm,'converged')=converged
		return(betaglm)
	}

	#Initialize some basic variables
	X0=X
	n=nrow(X)
	p=ncol(X)
	f=family()
	no.penalty=which(penalty.factor==0)
	if(length(no.penalty)==0) no.penalty=NULL
	
	#Check for warnings
	type=match.arg(type)
	if(type=='lars'&f$family!='gaussian'){
		warning("type lars is just for gaussian family; set to glmnet")
		type='glmnet'
	}
	if(type=='lars'&!is.null(no.penalty)){
		warning("type lars does not allow no penalty subsets; set to glmnet")
		type='glmnet'
	}
	if(type=='lars'&!is.null(offset)){
		warning("type lars does not allow offset; set to glmnet")
		type='glmnet'
	}
	
	#Perform one iteration to estimate beta0 when family is not Gaussian, and there is an intercept and/or there are unpenalized columns, and beta0 is unknown
	if(beta0[1]=='iter'&!(!intercept&is.null(no.penalty))&f$family!='gaussian'){
		#get qut estimates with beta0=glm(y~A)
		outqut=qut(y=y,X=X,fit=fit,family=family,alpha.level=alpha.level,M=M,qut.standardize=qut.standardize,
			intercept=intercept,offset=offset,bootstrap=bootstrap,sigma=sigma,beta0='noiter',estimator=estimator,
			type=type,lambda.seq=lambda.seq,penalty.factor=penalty.factor,lambda.min.ratio=lambda.min.ratio,nlambda=nlambda,lambda=lambda,...)	
		#Get new value for beta0
		if(intercept) indbeta0=c(1,no.penalty+1) else indbeta0=no.penalty+1
		beta0=outqut$beta[indbeta0]
	}
	else if(!is.numeric(beta0)) beta0=NA
	
    #Obtain lambdaqut (Quantile Universal Threshold)
	outqut=lambdaqut(y=y,X=X,alpha.level=alpha.level,M=M,qut.standardize=qut.standardize, family=family,intercept=intercept,no.penalty=no.penalty,offset=offset,bootstrap=bootstrap,beta0=beta0)
	X=outqut$Xnew #Get standardized matrix
	lambdaqut=outqut$lambda
	
	#Check if lambda is Inf: Problems estimating the model
	if(lambdaqut==Inf){
		beta=rep(0,n+1)
		betaglm=rep(0,n+1)
		beta[1]=NA
		betaglm[1]=NA
		warning('Difficulty estimating this model, try increasing M')
		if(type=='glmnet') lambdamax=outqut$lambda.max*sum(penalty.factor)/p
		else lambdamax=outqut$lambda.max*n
	}
	
	#No problems with lambda
	else{
		#Estimate sigma for gaussian family (if not specified)
		if(f$family=='gaussian'){
			if(!is.numeric(sigma)){
				if(sigma=='ols'){
					if(n<p){
						warning("since n<p, is not possible to calculate the OLS sigma estimator; set to 'qut'")
						sigma=sigmaqut(y=y,X=X,intercept=intercept,estimator=estimator,alpha.level=alpha.level,M=M,qut.standardize=qut.standardize,offset=offset,penalty.factor=penalty.factor,...)
					}
					else{
						if(intercept) sigma=sqrt(sum(lm(y~X,offset=offset)$res^2)/n)
						else sigma=sqrt(sum(lm(y~X-1,offset=offset)$res^2)/n)
						if(sigma<1e-9){
							warning("overfitting data for sigma estimation; set to 'qut'")
							sigma=sigmaqut(y=y,X=X,intercept=intercept,estimator=estimator,alpha.level=alpha.level,M=M,qut.standardize=qut.standardize,offset=offset,penalty.factor=penalty.factor,...)
						}
					}
				}
				else if(sigma=='qut') sigma=sigmaqut(y=y,X=X,intercept=intercept,estimator=estimator,alpha.level=alpha.level,M=M,qut.standardize=qut.standardize,offset=offset,penalty.factor=penalty.factor,...)
				else{
					warning("sigma rcv or cv does not take into account penalty factor")
					if(sigma=='rcv') sigma = sigmarcv(X=X,y=y,cv=FALSE,intercept=intercept)$sigmahat
					else if(sigma=='cv') sigma = sigmarcv(X=X,y=y,cv=TRUE,intercept=intercept)$sigmahat
					else stop('sigma must be a positive real number or a valid estimator type')
				}
			}
			else if(sigma<=0) stop('sigma must be a positive real number or a valid estimator type')
		}
		else sigma=1
		
		#Fit model
		if(type=='glmnet'){ #GLMNET
			lambdamax=outqut$lambda.max*sum(penalty.factor)/p
			#set sequence of lambdas
			if(lambda.seq==2){  #default glmnet values
				if(missing(lambda.min.ratio)) lambda.min.ratio=ifelse(n<p,0.01,0.0001)
				if(missing(nlambda)) nlambda=100
				if(missing(lambda)) lambda=NULL
			}
			else{ #starting from lambdaqut to lambdamax
				lambda.min.ratio=lambdaqut*sigma/lambdamax*0.9999

				if(lambda.seq==0 ){ #equi-spaced
					lambda=NULL
					if(lambda.min.ratio>1) lambda.min.ratio=0.9999 
				}
				else if(lambda.seq==1){ #log-scale
					if(lambda.min.ratio>1)	lambda=exp(seq(log(lambdamax),log(lambdamax*0.9999),length=nlambda))/n #n lambdas in logscale
					else lambda=exp(seq(log(lambdamax),log(lambdaqut*sigma),length=nlambda))/n #n lambdas in logscale
				}
			}
			if(missing(penalty.factor)) penalty.factor=rep(1,p)
			if(!missing(no.penalty)) penalty.factor[no.penalty]=0

			#Obtain GLM L1 path using GLMNET
			if(missing(fit)) fit=glmnet(X,y,standardize=FALSE,intercept=intercept,family=f$family,penalty.factor=penalty.factor,offset=offset,nlambda=nlambda,lambda.min.ratio=lambda.min.ratio,lambda=lambda,...)		

			if(max(fit$lambda)==Inf) beta=rep(0,p+1)
			else beta=coef(fit, s=lambdaqut*sigma/n*sum(penalty.factor)/p,offset=offset)/c(1,outqut$scale.factor) 

			lambdaqut=lambdaqut/n

			#Get GLM fitted coefficients
			betaglm=try(betaGLMestimates(),silent=TRUE)
		
		}
		else if(type=='lars'){ #LARS
		
			#LASSO fit
			lambdamax=outqut$lambda.max*n
			if(missing(penalty.factor)) penalty.factor=rep(1,p)
			X=X/penalty.factor
			
			if(missing(fit)) fit=lars(X,y,intercept=intercept,normalize=FALSE,...)
			beta=coef(fit,s=lambdaqut*sigma,mode='lambda')
			
			if(intercept) beta0=mean(y)-mean(X%*%beta) else beta0=0
			beta=beta/penalty.factor/outqut$scale.factor
			
			#Get GLM fitted coefficients
			betaglm=try(betaGLMestimates(),silent=TRUE)
			
			beta=c(beta0,beta)

			names(beta)[1]='(Intercept)'	
			names(beta)[-1]=1:p	
			
		}
	}

	#OUTPUT
	out=NULL
	out$lambda=lambdaqut
	out$lambda.max=lambdamax/n
	out$scale.factor=outqut$scale.factor
	out$beta=beta
	if(class(betaglm)[1]=="try-error") warning("No valid GLM fits were possible to find")
	else betaglm[is.na(betaglm)]=0
	out$betaglm=betaglm
	out$fit=fit
	out$family=family
	if(f$family=='gaussian') out$sigma=sigma
	else out$sigma=NA
	class(out)='qut'
	return(out)
}


