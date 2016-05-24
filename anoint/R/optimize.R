# CREATE DATA FRAME FOR EXACT PIM WITH THETA MULTIPLE APPLIED

data.frame.pim <- function(formula.anoint,data,theta=1){
	
	X <- model.matrix(update(formula.anoint@prognostic,
						paste("~.+",formula.anoint@trt,sep="",collapse="")),data)
						
	if(colnames(X)[1]=="(Intercept)") X <- X[,-1]

	trt.index <- which(colnames(X)==formula.anoint@trt)
	X[X[,trt.index]==1,-trt.index] <- X[X[,trt.index]==1,-trt.index]*theta
	
	df <- as.data.frame(X)
	
	response.index <- ncol(cbind(model.frame(formula.anoint@prognostic,data)[,1]))
	response <- subset(data,select=all.vars(formula.anoint@prognostic)[1:response.index])
	df <- cbind(df,response)

df
}


pim.exact.fit <- function(formula.anoint,data,theta){
	 # RETURNS RESTRICTED MLE FIT GIVEN PROPORTIONAL INTERACTION TERM
	fit.data <- data.frame.pim(formula.anoint,data,theta)
	
	if(formula.anoint@family=="coxph"){
		fit <- coxph(formula.anoint@prognostic.trt,data=fit.data)
	}
	else{
		fit <- glm(formula.anoint@prognostic.trt,data=fit.data,family=formula.anoint@family)
	}

fit
}

pim.exact.loglik <- function(theta,formula.anoint,data){
	 # RETURNS RESTRICTED MLE FIT GIVEN PROPORTIONAL INTERACTION TERM
	fit.data <- data.frame.pim(formula.anoint,data,theta)
		
	if(formula.anoint@family=="coxph"){
		fit <- coxph(formula.anoint@prognostic.trt,data=fit.data)
		-fit$loglik[2]
	}
	else{
		fit <- glm(formula.anoint@prognostic.trt,data=fit.data,family=formula.anoint@family)
		-as.numeric(logLik(fit))
	}
}

optimize.pim <- function(interval,formula.anoint,data,...){
	optimize(pim.exact.loglik,interval=interval,formula.anoint=formula.anoint,data=data,...)
}


pim.exact <- function(formula.anoint,data,interval=c(-5,5),...){
	
	result <- optimize.pim(interval=interval,formula.anoint,data,...)
	
	theta <- result$min
	fit <- pim.exact.fit(formula.anoint, data, theta=theta)

	alpha.index <- length(coef(fit))
	if(names(coef(fit))[1]=="(Intercept)") alpha.index <- c(1,alpha.index)
	
	loglik.null <- pim.exact.loglik(formula.anoint,data,theta=1)

	list(
		alpha = coef(fit)[alpha.index],
		vcov = vcov(fit),
		beta.control = cbind(coef(fit)[-alpha.index]),
		theta = theta,
		LRT = 2*(loglik.null-result$obj)
	)
}
