pim.fit <- function(formula, trt, data, family="binomial", na.action=na.omit,...){
	
		if(missing(trt))
			stop("Must indicate the treatment variable.")
			
		if(class(trt)!="character")
			stop("Treatment argument must be the character name of the treatment group variable.")

		data <- na.action(data[,unique(c(all.vars(formula), trt))]) # APPLY NA.ACTION
		
		# SEPARATE CONTROL AND TREATMENT DATASETS
		control <- data[data[,trt]==0,]
		treatment <- data[data[,trt]==1,]
		
		fitter <- function(formula, control, treatment, family){
			if(family=="coxph"){
				fit0 <- coxph(formula, data=control, ...)
				fit1 <- coxph(formula, data=treatment, ...)
				Sigma <- (vcov(fit0)+vcov(fit1))/2
				Sigma.root <- solve(chol(Sigma))
				U0 <- as.numeric(Sigma.root%*%fit0$coef)
				U1 <- as.numeric(Sigma.root%*%fit1$coef)
			}
		else{
				fit0 <- glm(formula, family=family, data=control, ...)
				fit1 <- glm(formula, family=family, data=treatment, ...)		
			if(length(grep("Intercept",fit0$coef))!=0){
				Sigma <- (vcov(fit0)[-1,-1]+vcov(fit1)[-1,-1])/2
				Sigma.root <- solve(chol(Sigma))
				U0 <- as.numeric(Sigma.root%*%fit0$coef[-1])
				U1 <- as.numeric(Sigma.root%*%fit1$coef[-1])
			}	
			else{
				Sigma <- (vcov(fit0)+vcov(fit1))/2
				Sigma.root <- solve(chol(Sigma))
				U0 <- as.numeric(Sigma.root%*%fit0$coef)
				U1 <- as.numeric(Sigma.root%*%fit1$coef)
				}
		}		
		
		fit <- propglm(U0, U1)
	
	list(
		interaction = as.numeric(fit$a),
		LRT = as.numeric(fit$T),
		lower= fit$lower,
		upper= fit$upper,
		pvalue = pchisq(as.numeric(fit$T),df=1,lower.tail=FALSE),
		model0 = fit0,
		model1 = fit1
	)
}
	
fitter(formula, control=control, treatment=treatment, family=family)
}

