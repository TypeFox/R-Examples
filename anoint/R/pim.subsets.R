anoint.onebyone.formula <- function(formula, trt){
	x <- as.character(formula)
as.formula(paste(x[2],x[1],"(",x[3],")","*",trt,sep="",collapse=""))
}

cor.propglm.value <- function(subset1, subset2){
	
	subset1 <- strsplit(subset1,"")[[1]]
	subset2 <- strsplit(subset2,"")[[1]]
	n1 <- sum(as.numeric(subset1))
	n2 <- sum(as.numeric(subset2))
	m <- sum(as.numeric(subset1)[subset1==subset2])
	
m/sqrt(n1*n2)
}

cor.propglm <- function(subsets){
	
	rho <- matrix(0,length(subsets), length(subsets))

	for(i in 1:nrow(rho)){
		for(j in 1:nrow(rho)){
			rho[i,j] <- cor.propglm.value(subsets[i], subsets[j])
		}
	}

mean(rho)
}


propglm <- function(U0,U1){
	
	R <- sqrt(t(U1)%*%U1)/sqrt(t(U0)%*%U0)
	Angle <- 	(t(U0)%*%U1)/(sqrt(t(U1)%*%U1)*sqrt(t(U0)%*%U0))
	Angle <- 2*Angle
	
	a <- (R-1/R)+sqrt((R-1/R)^2+Angle^2)
	a <- a/Angle
	
	u0.hat <- (a*U1+U0)/(a^2+1)
	u1.hat <- a*u0.hat
	        
     T0 <- t(U1-U0)%*%(U1-U0)/2
	 T1 <- t(U1-u1.hat)%*%(U1-u1.hat) # ALTERNATIVE
	 T2 <- t(U0-u0.hat)%*%(U0-u0.hat) # ALTERNATIVE
	 T <- T0-T1-T2

	roots <- function(a,b,c) {
		if(b^2<4*a*c)
			NA
		else
			(-b+c(1,-1)*sqrt(b^2-4*a*c))/(2*a)
	}
	
	a0 <- sum(U0^2)
	b0 <- -2*t(U0)%*%U1
	c0 <- sum(U1^2)
	k <- qchisq(.95,df=length(U0))
	
#	ci95.quad <- roots(a0-k,b0,c0-k)	# DEPRECATED
#	ci95.quad <- sort(ci95.quad)
	
	# BOOTSTRAP
	f <- function(x) ((x-1/x)+sqrt((x-1/x)^2+Angle^2))/Angle
	f <- Vectorize(f)
	p <- length(u0.hat)
	
	ncp0 <- t(u0.hat)%*%u0.hat 
	ncp1 <- t(u1.hat)%*%u1.hat 
	R <- sqrt(rchisq(n=1000,df=p,ncp=ncp1)/rchisq(n=1000,df=p,ncp=ncp0))
	ci95 <- quantile(f(R),c(.025, .975))
	
list(a=a,T=T, lower=ci95[1], upper=ci95[2])
}

getU <- function(model){
	
	if(length(grep("Intercept",names(model$coef)))==1){ # REMOVE INERCEPT
		Sigma <- chol(vcov(model)[-1,-1])
		as.numeric(solve(Sigma)%*%model$coef[-1])
		}
	else{
		Sigma <- chol(vcov(model))
		as.numeric(solve(Sigma)%*%model$coef)
	}
}


subsets.index <- function(p){
	
	# p = number of covariates
	n.subsets <- 2^p
	subset.index <- matrix(0,ncol=p,nrow=n.subsets)
	
	repeats <- 2^(0:(p-1))
	repeats <- repeats[p:1]
	
	for(i in 1:p){
		subset.index[,i] <- rep(rep(c(0,1),each=repeats[i]),lengt=n.subsets)
	}

subset.index
}
	
propglm.subsets <- function(U0,U1,alpha=0.05){
	
	subset.index <- subsets.index(length(U1))[-1,] # EXCLUDE NULL MODEL
	
	subset.index.string <- apply(subset.index,1,function(x) paste(x,sep="",collapse=""))
	
	which.index <- apply(subset.index,1,function(x) which(x==1)) # LIST OF ALL SUBSETS
	
	LRT <- sapply(which.index, function(index) propglm(U0[index],U1[index]))
	
	result <-  data.frame(
		subset = subset.index.string,
		interaction = as.numeric(LRT[1,]),
		LRT = as.numeric(LRT[2,]),
		lower= as.numeric(LRT[3,]),
		upper= as.numeric(LRT[4,]),
		pvalue = pchisq(as.numeric(LRT[2,]),df=1,lower.tail=FALSE)
	)

	result <- result[order(result$LRT,decreasing=TRUE),]	 # SORT
	result$subset <- as.character(levels(result$subset))[result$subset]
	rho.bar <- cor.propglm(result$subset)
	result$reject <- result$pvalue<=alpha/((1-rho.bar)*nrow(result))

result
}



pim.subsets <- function(formula, trt, data, family="binomial", na.action=na.omit, fwer=0.05,...){
	
		if(missing(trt))
			stop("Must indicate the treatment variable.")
			
		if(class(trt)!="character")
			stop("Treatment argument must be the character name of the treatment group variable.")
			
		data <- na.action(data[,unique(c(all.vars(formula), trt))]) # APPLY NA.ACTION
				
		# SEPARATE CONTROL AND TREATMENT DATASETS
		control <- data[data[,trt]==0,]
		treatment <- data[data[,trt]==1,]
		
		# FIT MODELS
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

propglm.subsets(U0=U0, U1=U1, alpha=fwer)
}


anoint.subgroups <- function(formula, trt, data, family="binomial", na.action=na.omit, fwer=0.05,...){

		if(missing(trt))
			stop("Must indicate the treatment variable.")
			
		if(class(trt)!="character")
			stop("Treatment argument must be the character name of the treatment group variable.")
			
		covariates <- all.vars(update(formula,NULL~.))
		data <- na.action(data[,unique(c(all.vars(formula), trt))]) # APPLY NA.ACTION
		
		# SUBGROUP ANALYSIS
		anoint.fit <- anoint(anoint.onebyone.formula(formula, trt), data=data, family=family)
		onebyone.results <- obo(anoint.fit)
		subgroup.lrts <- onebyone.results$LRT
		subgroup.pvalues <- onebyone.results$pvalue
		n <- ifelse(family=="coxph",3,4)
		cov.n <- n-2
		q <- qnorm(1-fwer/2)

		subgroup.ints <- sapply(onebyone.results$fit, function(x) coef(x)[n])
		subgroup.width <- sapply(onebyone.results$fit, function(x) sqrt(vcov(x)[n,n])*q)
		subgroup.lower <- subgroup.ints-subgroup.width
		subgroup.upper <- subgroup.ints+subgroup.width
	    subgroup.index <- apply(diag(length(covariates)),1, function(x) paste(x,sep="",collapse=""))
		include.exclude.matrix <- diag(length(covariates))==1
		
		result <- data.frame(
		    subset = subgroup.index,
			interaction = subgroup.ints,
			LRT = subgroup.lrts,
			lower= subgroup.lower,
			upper= subgroup.upper,
			pvalue = subgroup.pvalues,
			covariates = covariates,
			reject = subgroup.pvalues<=fwer/length(covariates)
	)
	
	row.names(result) <- covariates
	o <- order(subgroup.lrts, decreasing=TRUE)	# SORT BY LARGEST LRT
	result <- result[o,]

	print(result)

	result <- as.list(result)
	result$include.exclude.matrix <- include.exclude.matrix[o,]
	
	invisible(result)
}



pim.subsets <- function(formula, trt, data, family="binomial", na.action=na.omit, fwer=0.05,...){
	
		if(missing(trt))
			stop("Must indicate the treatment variable.")
			
		if(class(trt)!="character")
			stop("Treatment argument must be the character name of the treatment group variable.")
			
		covariates <- all.vars(update(formula,NULL~.))
		data <- na.action(data[,unique(c(all.vars(formula), trt))]) # APPLY NA.ACTION
		
		# SEPARATE CONTROL AND TREATMENT DATASETS
		control <- data[data[,trt]==0,]
		treatment <- data[data[,trt]==1,]
		
		# FITTER FUNCTION
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
			
				fit0 <- glm(formula, family=family, data=control)
				fit1 <- glm(formula, family=family, data=treatment)		
		
		if(length(grep("Intercept",names(fit0$coef)))!=0){ # INTERCEPT CHECK
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
		
		propglm(U0=U0, U1=U1)
}
	
	subset.index <- subsets.index(length(covariates))[-1,] # EXCLUDE NULL AND SUBGROUP MODELS    

    subgroup.index <- rowSums(subset.index)==1
	subset.index <- subset.index[!subgroup.index,]

	if(!is.matrix(subset.index)) subset.index <- matrix(subset.index, nrow=1)

	subset.index.string <- apply(subset.index,1, function(x) paste(x,sep="",collapse=""))
	include.exclude.matrix <- subset.index==1
	
	# CREATE FORMULA LIST
	formula.list <- apply(subset.index,1, 
		function(x) update(formula,paste("~",paste(covariates[x==1],collapse="+")),collapse=""))

	fits <- sapply(formula.list, function(f) fitter(f,control=control,treatment=treatment,family=family))
	
	result <-  data.frame(
		subset = subset.index.string,
		interaction = as.numeric(fits[1,]),
		num = apply(subset.index,1,sum),
		LRT = as.numeric(fits[2,]),
		lower= as.numeric(fits[3,]),
		upper= as.numeric(fits[4,]),
		pvalue = pchisq(as.numeric(fits[2,]),df=1,lower.tail=FALSE)
	)
	
	result$subset <- as.character(levels(result$subset))[result$subset]
	rho.bar <- cor.propglm(result$subset)
	result$reject <- result$pvalue<=fwer/((1-rho.bar)*nrow(result))
	o <- order(1-result$reject, result$LRT,decreasing=FALSE) # SELECT LOWEST LRT OF SELECTED
	result <- result[o,]	 # SORT
	
	print(result)
	
	result <- list(
	    subset = result$subset,
		interaction = result$interaction,
		LRT = result$LRT,
		num = result$num,
		lower= result$lower,
		upper= result$upper,
		pvalue = result$pvalue,
		include.exclude.matrix = include.exclude.matrix[o,],
		covariates = covariates,
		reject = result$reject
	)

	invisible(result)
}
