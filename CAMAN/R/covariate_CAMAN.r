#
#`!`<- function(x)
#	if (inherits(x, "character") == FALSE)
#				.Primitive("!")(x) else invisible(x)

mixcov <- function(dep,fixed,random="",data,k,weight=NULL, pop.at.risk=NULL, 
                   var.lnOR=NULL, family="gaussian", maxiter=50, 
                   acc=10^-7, returnHomogeneousModel = FALSE)
{
  stopifnot(family %in% c("gaussian", "poisson")) # disable binomial for now as it is not working as intended
    cl <- match.call()
	nn <- nrow(data)	
	if (!is.null(pop.at.risk) && family=="binomial"){
		binomDat <- cbind(observed=data[[which(colnames(data)==dep)]], 
                      popRisk=data[[which(colnames(data)==pop.at.risk)]]-data[[which(colnames(data)==dep)]])
		dep="binomDat"
	}
	
	if (is.null(pop.at.risk)) pop.at.risk <- rep(1, nn)
	else {
		if (is.character(pop.at.risk) ) pop.at.risk = data[,which(colnames(data)==pop.at.risk)]
	}
	if (is.null(var.lnOR)) var.lnOR <- rep(1, nn)
	else {
		if (is.character(var.lnOR) ) var.lnOR = data[,which(colnames(data)==var.lnOR)]
	}
	#if (is.null(weight)) pop.at.risk <- rep(1, nn)


    # Usual homogenous model 
    form<-as.formula(paste(paste(dep,"~"),paste(fixed,collapse="+"))) #dependencies for linear model
    m<-glm(form,family=family,weights=weight,data=data,x=T,na.action=na.omit, offset=log(pop.at.risk)) #compute linear model

    #y0 <- as.vector(model.extract(model.frame(formula(m), data = data), response))
	ixDep <- which(names(data) == dep)
	ixFixed <- which(names(data) %in% fixed)
	ixRandom <- which(names(data) %in% random)
	y <- data[,ixDep]
	y0 <- y
	ixColSort <- c(ixDep, ixFixed, ixRandom )
    
	idxControl_dep <- 1; names(idxControl_dep) = dep;
	idxControl_fixed <- integer(0)
	idxControl_Random <- integer(0)
	
	if(length(ixFixed)>0) {idxControl_fixed <- 2:(2+length(ixFixed)-1); names(idxControl_fixed) = names(data)[ixFixed];}	
	if(length(ixRandom)>0) {idxControl_Random <- 2:(2+length(ixRandom)-1); names(idxControl_Random) = names(data)[ixRandom];}
	
	idxControl <- list(ixDep = idxControl_dep, ixFixed=idxControl_fixed, ixRandom=idxControl_Random)
	colnmes = names(data)[ixColSort]
	data <- as.data.frame(data[,ixColSort])  #rearrange data sothat the dependent variable is in the first column
	#idxControl has the indices of the rearranged data!!
	names(data) = colnmes
	mix_data <- data.frame(data[,1], rep(1,nn), pop.at.risk, var.lnOR)
	#no weights
	
    if(is.null(weight) && family=="gaussian")
    {
        dfg1 <- m$df.residual
        wt0<- deviance(m)/dfg1
    }
    else wt0<-1./weight

    pre<-fitted(m) #generate data with homogeneous model m
    
    #compute Log-Likelihood of the homogeneous model 
	logl0=0
    if (family=="gaussian")
        logl0<-sum(log(dnorm(y0,pre,sqrt(wt0))))
    else if (family =="poisson")
        logl0<-sum(log(dpois(y0,pre)))
        
#    cat("\n", "Fit of the usual model:", "\n")
#            print(summary(m))
#            cat("\n", "resid.variance:", round(wt0, 5))

	BIC_homo <- -2*logl0+log(nn)*length(coef(m))
	LL_homo <-   logl0
    
    #print(logl0)
    if(k==1) 
    {
    return(list(m, BIC=BIC_homo, LL=LL_homo)) #only one component
    }

	x0<-m$x
    #form<-as.formula(paste(paste(dep,"~"),paste("Z")))
    form <- paste(paste(dep,"~",sep=""),paste("Z+",collapse="+"))
    form <- paste(form,paste(fixed,collapse="+")) 
    ### form <- paste(fixed[fixed!="1"], collapse="+")
    if((random[1] !="") || (length(random)>1)){
		random <- random[random!="1"] #cut out the intercept 
        form_rand<-paste("Z/",random,sep="",collapse="+")
        form<-paste(form,form_rand,sep="+")
    }
    form<-as.formula(paste(form,"-1"))
 
    # starting values
    p <- rep(1/k,k) #components are equal distributed 
	ixCols_fixedEffects = which(colnames(data) %in% fixed);
	ixCols_randomEffects = which(colnames(data) %in% random);
	
    if ((length(fixed)==1)&&(fixed =="1")) numPara_fixed<-0  #no fixed effects
	else if (sum(sapply(data[,ixCols_fixedEffects], is.factor))== 0) numPara_fixed<-length(fixed) ##?? numPara_fixed<-length(fixed)-1 or sum(fixed!="1") 
	else{
		if (length(ixCols_fixedEffects)==1) ixFactorsFixed = which(lapply(data, is.factor)[[ixCols_fixedEffects]])
		else ixFactorsFixed = which(unlist(lapply(data, is.factor)[ixCols_fixedEffects]))
		#if clause is just needed to seperate between [[single_idx] and [several_idx]
		#number_of_levels_Fixed = unlist(lapply(data, function(xx){length(levels(xx))}))[ixCols_randomEffects[ixFactorsRandom]]
		number_of_levels_Fixed = unlist(lapply(data, function(xx){length(levels(xx))}))[ixCols_fixedEffects[ixFactorsFixed]]  
		numPara_fixed <- (length(fixed)-length(ixFactorsFixed)) + sum(number_of_levels_Fixed -1)
	}
	if((length(random) == 1) && (random == "")) numPara_random <- 0 #no random effects
    else if (sum(sapply(data[,ixCols_randomEffects], is.factor))== 0) numPara_random<-k*length(random) 
			#--> there were just numeric (no factorized) data, so, numPara_random = k*nRandom
	else{ if (length(ixCols_randomEffects)==1) ixFactorsRandom = which(lapply(data, is.factor)[[ixCols_randomEffects]])
		else ixFactorsRandom = which(unlist(lapply(data, is.factor)[ixCols_randomEffects]))
		#if clause is just needed to seperate between [[single_idx] and [several_idx]
#		cat("ixFactorsRandom", ixFactorsRandom)
		number_of_levels_Random = unlist(lapply(data, function(xx){length(levels(xx))}))[ixCols_randomEffects[ixFactorsRandom]]
		numPara_random <- k*(length(random) - length(ixFactorsRandom)) #no factors ... just k parameter for a factor
		numPara_random <- numPara_random + sum((number_of_levels_Random-1)*k)
	 # -1 needed because in the resulting Model, the last factor is considered to be the intercept
	 }
    numPara<-numPara_fixed+k+numPara_random  #total no. of parameters

    b<-rep(0,numPara) #parameter initialization
	b <- seq(min(y0), max(y0), length.out=k)
	## b[1:k] -> estimated intercepts of the components (starting values)
	if (numPara>k) b[(k+1):numPara] = 0
	## b[k+1:numPara] -> starting values for other parameters
	
    # obtain solution of EM-algorithm
    mem<-mix.perform_glm(form,data,k,p,y,b,var=wt0,family=family, 
                         maxiter=maxiter, acc=acc, expected = pop.at.risk)
    logem<-mem$logl

    p <- mem$p
    #wp <- p[iii]
    x <- mem$x
    xf <- mem$xf
    m1 <- mem$m1 
	steps <- mem$n_iter

#    cat("\n", " Fit of the ", round(k, 2), "-component mixture model:", "\n")
#    cat("\n", "coefficients:", "\n")
	coefMatrix <- coefficients(summary(mem$m1))
	coefMatrix[, 2] <- coefMatrix[, 2] * sqrt(2)
	coefMatrix[, 3] <- coefMatrix[, 1]/coefMatrix[, 2]
    
	#extract information out of the coefMatrix. We make use of the property 
	#that the first rows are the intercepts, then the fixed effects are  
	
	if (nrow(coefMatrix)!= numPara)
		warning("NUMBER OF PARAMETERS SEEMS TO BE INCONSISTENT!!")
	
	is.intercept = rep(NA, nrow(coefMatrix)); is.fixed = rep(NA, nrow(coefMatrix)); 
	is.random = rep(NA, nrow(coefMatrix));
	#if the value is unequal to NA, the value indicates the membership to a component 
	is.intercept[1:k] = 1:k
	if (numPara_fixed>0) is.fixed[(k+1):(k+numPara_fixed)] = TRUE  # member of all components! (--> fixed)
	if (numPara_random >0 ) 
		is.random[(k+numPara_fixed+1):(k+numPara_fixed+numPara_random)] = rep(1:k, length.out=numPara_random)
	coefMatrix <- data.frame(coefMatrix, is.intercept=is.intercept, is.fixed=is.fixed, is.random=is.random)
	coef<-coef(mem$m1)
	meancoef <- 0
	for (i in 1:k) 
		meancoef <- suppressWarnings(meancoef + p[i] * (coefMatrix[which(coefMatrix$is.intercept==k),1] + 
		    sum(coefMatrix[which(coefMatrix$is.random==k),1] * as.numeric(mean(data[idxControl$ixRandom])), na.rm=TRUE) + sum(coefMatrix[which(coefMatrix$is.fixed),1] * as.numeric(mean(data[idxControl$ixFixed])) , na.rm=TRUE)))
	#we need to suppress the warnings due to factorial data. mean(factors) == NA, 
	#thus here, factorial data doesn't have any influence in the common effect!

	hetvar <- 0
	for (i in 1:k)
		hetvar = hetvar + p[i]*((meancoef - (coefMatrix[which(coefMatrix$is.intercept==k),1] + 
		   mean(coefMatrix[which(coefMatrix$is.random==k),1]) + mean(coefMatrix[which(coefMatrix$"is.fixed"),1])))^2)
	
	hetvar <- sqrt(mean(rep(p,nn) * (as.vector(fitted(m1)) - meancoef)^2))	# TODO check wheter hetvar correct.. sqrt seems to be false
	#meancoef <- mean(as.vector(fitted(m1))[classification_tmp])    # TODO check wheter meancoef correct
	
	residVar <- as.numeric(NA)
	if (family=="gaussian"){ 
		residVar <- mem$residVar
	}
	#hetvar = sum(p * (meancoef - coef)^2)
	
	
	
	pPosteriori<-matrix(mem$pPosteriori,nrow=nn,ncol=k)
	
	resultObj <- new("CAMAN.glm.object",dat=mix_data, family=family, LL=mem$logl, 
			num.k=k, p=p, t=coef[1:k], num.obs = nn, steps=steps, 
			otherParams = c(maxiter, acc), BIC=-2*logem+log(nn)*(length(b)+k-1), 
			commonEffect = meancoef, hetvar = hetvar, coefMatrix=coefMatrix, 
			numPara=numPara, cl= cl, depVar=dep, fixedVar = fixed, random = random,
			form=form, glmModel=m1, residVar = residVar, idxControl=idxControl, inputData=data)
	posterior_matrix <- mix.densDistr(resultObj)
	resultObj@classification = as.numeric(apply(posterior_matrix, 1, which.max)) #as.numeric(apply(posterior_matrix, 1, which.max))
	resultObj@prob <- posterior_matrix# posterior_matrix
	resultObj@fittedObs <- fitted(m1)[resultObj@classification + ((0:(nn-1))*k)]
#	cat("fittedObs", resultObj@fittedObs)
	
	#resultObj@hetvar = (0:(nn-1))*k + resultObj@classification
	
#	cat("sumS1",sum(mem$s1))
    
	if (returnHomogeneousModel)
		return(list(mixModel = resultObj, homoModel = list(lm=m, LL=LL_homo, BIC=BIC_homo) ) )
    return(resultObj)
}



########################################

mix.perform_glm <- function(form,data,k,p=NULL,y,b=NULL,
                            expected=NULL, var=NULL,weight=NULL,family="gaussian", 
                            shuffle=FALSE, maxiter=30, acc=10^-7)
  
## PD: argument shuffle is redundant right now, but doesn't hurt
{
	nn<-length(data[,1])  #length of !!non-expanded!! data
# data augmentation
	ii<-rep(1:nn,each=k)#e.g.:k = 3 --> ii = 111222333444...
	#y<-y0[ii]
	iii<-rep(1:k,nn)#e.g.:k = 3 --> ii = 12341234...

	
	# model
	#i<-rep(1:nn,each=k)#same as ii: e.g.:k = 3 --> ii = 111222333444...
	dataExpanded <- data.frame(data[ii,])		 #expand the data with factor k
	names(dataExpanded) <- names(data)
	expected = expected[ii]
	dataExpanded$Z<-as.factor(rep(1:k,nn)) #add a categorial variable for group-assignment
	y<-y[ii] #expand the outcome
	#set weights
	#if(is.null(weight)) wt<-rep(1,nn)
	#else wt<-residVar[ii]
		
	#set initial weights
	if(is.null(weight)){
	       residVar<-var
	}
	else residVar<-1/weight
	
	
	ii<-rep(1:nn,each=k)  #get the indices 111222333
	wt<-rep(1,nn)    #set weigths
	
	# create starting values
	grad<-rep(10,k)
	
	iii<-rep(1:k,nn) #123412341234....
	
	pPosteriori  <- p[iii] #expand component weigths
	pre <- b[iii] 
	delta_acc <- 10
	continue_iterate <- TRUE
	n_shuffled <- 0
	#while (delta_acc > 10^-6){
	while (continue_iterate){
		n_shuffled <- n_shuffled +1
		#compute matrix of densities of the corresponding distribution
		xf <- computeDensities(family, y, pre, residVar, k, nn)

		#calculate mixture density 
		s1<- apply(xf*pPosteriori, 2,sum) #mischverteilungsdichte%
		
		# calculate new weights
		#pPosteriori <- as.vector(pPosteriori*xf/s1[ii]) #posteriori
		pPosteriori <- as.vector(exp(log(pPosteriori)+log(xf)-log(s1[ii]))) #posteriori
		
		#wt --> studienspezifische Gewichte, die mit den daten ?bergeben werden (weight) 
		dataExpanded$wgt <- pPosteriori/wt
		dataExpanded$expect <- log(expected)
		wgt <- pPosteriori/wt ## avoids CRAN policy problems with "no visible binding for global variable"
		expect <- log(expected) ## dito
    
		# model
		if(family=="poisson") b[1:k]=log(b[1:k])
		
		diffLL<- 10 #starting value...
		n_iter <- 0 
		while (diffLL > acc)#
			{
			n_iter <- n_iter + 1
#			cat ("diffLL=", diffLL, "\n")
			# new coefficients

			m1<-glm(form,family=family,weights=wgt,
              start=b,x=TRUE,
              data=dataExpanded, offset=expect)
			b<-coef(m1)
			if(is.null(weight) && family=="gaussian"){
			  	residVar <- deviance(m1)/nn #residualvariance
			}
		
			pre<-fitted(m1)
			
			#compute matrix of densities of the corresponding distribution
			xf <- computeDensities(family, y, pre, residVar, k, nn)
			
			# new mixing weights	
			fitW <- glm(pPosteriori~Z, family = gaussian, data = dataExpanded, 
                  na.action = na.omit)
			wp <- predict(fitW)
			
			#calculate mixture density
			s1 <- apply(xf*wp, 2,sum)
			p <- wp[1:k]
			
			# calculate new weights
			pPosteriori <- as.vector(wp*xf/s1[ii])
			for(i in 1:k) grad[i] <- sum(xf[i,]/s1)/nn
			dataExpanded$wgt <- pPosteriori/wt
			logl <- sum(log(s1))
			diffLL <- abs(max(grad)-1)
#			cat("diffLL_after =", diffLL)
			logl0 <- logl
		}
		x <- m1$x
		# emgb
## PD: I removed the complete shuffle block since it isn't working.
# 		if (shuffle){
# 			if(is.null(weight))wt<-residVar
# 			
# 			s1<-apply(wp*xf,2,sum)
# 			tmax<-my_grad(y0,s1,wt,x0,kk=40,family)
# 			if(family=="gaussian")
# 				temp<-dnorm(y0,tmax,sqrt(wt))
# 			else if (family=="poisson")
# 				temp<-dpois(y0,tmax) 
# #			cat("\n","SA max", tmax,"\n")
# 			
# 			for(i in 1:k){
# 				llnu[i]=sum(log(s1-p[i]*xf[i,]+p[i]*temp))
# 			}
# 			ix<-which.max(llnu)
# 			coef[ix]<-tmax
# 			b<-coef(m1)
# 			pre<-x%*%b
# 			xf <- computeDensities(family, y, pre, residVar, k, nn)
# 			s1<- apply(xf*wp, 2,sum)
# 			tt<-exp(log(wp)+log(xf))
# 			logs1<-apply(tt,2,lse)
# 			
# 			wt<-residVar ###???
# 			residVar<-wt
# 			delta_acc <-  abs(oldlog - logl)
# 			if ((delta_acc < 10^-6) || (n_shuffled >= maxiter)) continue_iterate = FALSE  #convergence criterio fulfilled
# 			else continue_iterate = TRUE #convergence criterion not fulfilled --> continue!
# 			oldlog<-log1 
# 			}
# 		else continue_iterate = FALSE  #there was no shuffeling, so just make one iteration and 
    continue_iterate = FALSE ## PD: need to drop this line if shuffling is reimplemented!
	}
	#cat("residualvarianz: ", residVar)
suppressWarnings(return(list(m1 = m1,p=p,pPosteriori = pPosteriori,xf = xf,x=x,logl = logl,
                             residVar = residVar,n_iter = n_iter, s1 =s1)))
}




my_grad <- function(y,s1,wt,x,kk=20,family)
{
	nn<-length(y)
	
	min<-min(y)
	max<-max(y)
	if(family=="poisson")
	{
		max<-log(max)
		min<-log(min)
	}
#$step<-(max-min)/(kk-1)
#$t<-seq(min,max,by=step)
	grad<-matrix(0,kk,kk)
	t<-matrix(0,kk,kk)
	count<-0
	co<-rep(0,2)
	a<-seq(min,max,len=kk)
	b<-seq(0.05,0.1,len=kk)
	maxi<--1000
	for(j in 1:kk)
	{
		for(i in 1:kk)
		{
			co[1]<-a[i]
			co[2]<-b[j]
			
			eta<-x%*%co
			if(family=="gaussian")
				xf<-dnorm(y,eta,sqrt(wt))
			else if(family=="poisson")
				xf<-dpois(y,exp(eta))
			grad[j,i]<-sum(xf/s1)/nn
			if(grad[j,i] > maxi)
			{
				maxi<-grad[j,i]
				imax<-i
				jmax<-j
			}
		}
	}
	op <- par(bg = "white")
	persp(a,b,grad,theta=30,phi=40,col="lightblue")
	co[1]<-a[imax]
	co[2]<-b[jmax]
	return(co)
}

computeDensities <- function(family, y, pre, residVar, k, nn){
	if(family=="gaussian")
		xf <- matrix(dnorm(y, pre, sqrt(residVar)), nrow = k, ncol = nn)  #TODO sqrt(residVar)/k ??? 
	if(family=="poisson")
		xf <- matrix(dpois(y, pre), nrow = k, ncol = nn)
	return(xf)
}

# warum wird die vorletzte Iteration zur?ckgegeben?


# objectShow <- function(){
# 	cat("\n", "Mixing weights: ", round(p, 4), "\n")
# 	cat("\n", "Coefficients: ", round(p, 4), "\n")
# 	cat("\n", "Common effect: ", meancoef, "\n")
# 	cat("\n", "Heterogeneity variance: ", hetvar)
# 	cat("\n", "Heterogeneity STD: ", sqrt(hetvar))
# 	cat("\n", "log-likelihood at iterate:",logem,"\n")
# 	cat("\n", "BIC ",-2*logem+log(nn)*(length(b)+k-1),"\n")		
# }



mix.effectsOfCovariate <- function(obj, effectnme){ #for fixed effects
	ixdat <- obj@idxControl$ixFixed[which(names(obj@idxControl$ixFixed) == effectnme)]
	
	if (is.numeric(obj@inputData[,ixdat])){ #no 
		ParaOfFixedEffects <- obj@coefMatrix[which(rownames(obj@coefMatrix)==effectnme),1]
		res <- ParaOfFixedEffects * obj@inputData[,ixdat]
	}
	else{
		ixCoef <- which(substr(rownames(obj@coefMatrix),1,nchar(effectnme))== effectnme)
		res <- rep(0, obj@num.obs)
		for (i in ixCoef){
			tmp_factor <- substr(rownames(obj@coefMatrix)[i], nchar(effectnme)+1, nchar(rownames(obj@coefMatrix)[i]))
			ixWithFactor <- which(obj@inputData[,ixdat] == tmp_factor)
			res[ixWithFactor] <- obj@coefMatrix[i,1]
		}
	return(res)
	}
}

mix.densDistr <- function(mix){
	# computes the probability for each observation (1..n - row of mix@dat) belonging to each component (1..k)
	# returns a matrix of dimension n x k
	dat <- mix@dat[,1]	
	res <- matrix(ncol=mix@num.k, nrow=length(dat))
	p <- mix@p
	
	if (class(mix)[1] == "CAMAN.object"){
			
		if (mix@family == "gaussian") {
			mu <- mix@t
			mix.sd <- sqrt(mix@component.var)
			for (i in 1:mix@num.k) res[,i] <- sapply(dat,
						function(x){p[i]*dnorm(x, mu[i], mix.sd ) / sum(p*dnorm(x, mu,
											mix.sd ))})
		}
		if (mix@family == "binomial") {
			prob <- mix@t
			popAtRisk <- mix@dat[,3]
			for (i in 1:mix@num.k) res[,i] <- apply(cbind(dat, popAtRisk), 1,
						function(x){p[i]*dbinom(x[1], x[2], prob[i]) / sum(p*dbinom(x[1],
											x[2], prob))})
		}
		if (mix@family == "poisson") {
			lambda <-  mix@t
			popAtRisk <- mix@dat[,3]
			for (i in 1:mix@num.k) res[,i] <- apply(cbind(dat, popAtRisk), 1, 
						function(x){p[i]*dpois(x[1], x[2]* lambda[i]) / 
									sum(p*dpois(x[1], x[2] * lambda))})
		}
	}
	else if (class(mix)[1]== "CAMAN.glm.object"){
	obs.hat <- matrix(as.vector(fitted(mix@glmModel)), ncol=mix@num.k, byrow=TRUE)
		for (i in 1:mix@num.k){
			if (mix@family == "gaussian") {
				mu <- mix@coefMatrix[1:mix@num.k,1]
				mix.sd <- sqrt(mix@residVar)
				for (j in 1:mix@num.k) res[,j] <- apply(cbind(dat,obs.hat), 1, 
							function(x){p[j]*dnorm(x[1], x[j+1], mix.sd ) / sum(p*dnorm(x[1], x[-1],mix.sd ))})				
			}
			if (mix@family == "poisson") {
				# TODO: poisson mix.densdistr 
				lambda <-  mix@coefMatrix[1:mix@num.k,1]
				popAtRisk <- mix@dat[,3]
				for (j in 1:mix@num.k) res[,j] <- apply(cbind(dat, popAtRisk, obs.hat), 1, 
							function(x){p[j]*dpois(x[1], x[2]* x[j+2]) / 
										sum(p*dpois(x[1], x[2] * x[-c(1:2)]))})
			}
			if (mix@family == "binomial") {
				prob <- mix@t
				popAtRisk <- mix@dat[,3]
				for (i in 1:mix@num.k) res[,i] <- apply(cbind(dat, popAtRisk), 1,
							function(x){p[i]*dbinom(x[1], x[2], prob[i]) / sum(p*dbinom(x[1],
												x[2], prob))})
			}
		}

	} 
	return(res)
	}
	
	
	
	summary.CAMAN.glm.object <- function(object, ...){
		cat("Summary of a Computer Assisted Mixture Analysis with covariates: \n \n")
		cat("Data consists of", object@num.obs, "observations (rows). \n")
		cat("The Mixture Analysis identified", object@num.k, "component")
		if (length(object@num.k) >0) cat("s")
		cat(" of a", object@family, "distribution: \n \n")
		n <- object@num.obs    
		cat("mixing weights:\n")
		p_tmp <- object@p
		names(p_tmp) = paste("comp.", 1:object@num.k)
		print(p_tmp)
		
		cat("\n Coefficients :\n")
		coefPrint <- object@coefMatrix[,1:4]
		print(coefPrint)
		if (object@family=="gaussian") cat("residual variance:", object@residVar)
		cat("\n")
		
		cat("\n Log-Likelihood:",object@LL,"    ")
		cat("BIC:",object@BIC,"\n")
		
		cat("regression formular: ")
		print(object@form)
		
		cat("iteration steps:", object@steps,"\n")
		cat("heterogeneity variance:", object@hetvar,"\n")
		cat("common effect:", object@commonEffect,"\n \n")
		
		cat("overview over the variables of the model:\n")
		cat("dependent variable:", object@depVar,"\n")
		cat("fixed effects:", object@fixedVar,"\n")
		cat("random effects:", object@random,"\n")
		
	}
		