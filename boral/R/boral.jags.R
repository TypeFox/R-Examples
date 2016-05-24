##############
## Latent Variable Ordination and Regression using MCMC 
## Site effects fitted as fixed effects
## NB parameterized as V = mu + phi*mu^2
## Ordinal data handled as propotional odds regression, with same cutoff points for all spp. 
## Multinomial regression is available, but CURRENTLY NOT ALLOWED DUE TO COMPUTATION TIME
## This should make sense if the same ordinal scale is applied to all species

## Changes from v0.7
## 1) Corrected MAJOR error for normal distribution random row effect!!!
## 2) par resets for lvsplot and plot.boral
## 3) Compound Laplace Metropolis ICs have been removed due to its instability
## 5) Probit link now used with binomial and ordinal family
## 6) Overdispersion parameters not sampled when they don't need to be, i.e. if all columns of y are binomial, exponential, multinomial, ordinal....
## 7) For reasons of computation and increased customizability, make.jagsboralmodel and make.jagsboralnullmodel now incoporate an argument for trial.size, specified in the same manner as in boral.default
## 8) make a note about sampling only one chain, the sign switching problem, and how it can screw up your ordination
## 9) Negative binomial parameterized in terms of size rather than dispersion
## 4) If data is entirelly Bernoulli, then a "bernoulli" family is used with the step parameterization and probit link, and the latent variables are constructed so that the variance of lambda%*%t(lambda) + psi = 1
## 10) Spp-specific coefficients and intercepts can now be regressed against traits, beta_0j and beta_jk ~ N(gamma_k*traits_j, sigma2_k). Create life also allows simulating responses from traits

## Changes from v0.8
## 1) qqplot in plot.boral now adds a 1:1 line for better diagnosis
## 3) get.enviro.cor and ger.residual.cor now return covariance and correlation matrices
## 4) get.residual.cor and ger.enviro.cor now also print out correlation matrices that contains only significant correlations (HPDintervals not containing zero). As a result, there is a prob argument now included in both, with a default of 0.95
## 5) Default settings for n.burnin, n.iteration, and n.thin have been ramped up
## 6) In lvsplot, the new.plot argument has been removed and replaced with a est argument to plot either posterior mean or median. Also, a main argument has been included to allow custom titles
## 7) For lvsplot, include an alpha parameter that allows the user to tweak the scaling themselves. This is akin to the different type of scalings in CA. For biplots, we generally prefer adjust alpha until the LVs and coefs are the same scale. This is usually around the value of alpha = 0.5 to 0.55
## 8) After some complaints from users not the swtich to size for the negative.binomial family, it has now been switched back to V = mu + phi*mu^2. However,...
## 9) All random effects distributions are now parametrized in terms of their standard deviation, with a uniform prior now placed on the standard deviation sigma2 instead of on the variance. This is following advice by Gelman, 2006, Prior distributions for variance parameters in hierarchical models
## 10) Similarly, all normal and lognormal responses are now parameterized in terms of standard deviation instead of variance.


## TODO: 1) correct problems with increasing number of LVs causing increase in covariance and variance i.e., scale problem?; 2) Calculate proprtion of deviance explained based on marginal or conditional log-likelihood? 3) Allow options for weakly informative priors in terms t/Cauchy and half-t distributions; 4) draw coefficients as random effects, or reduce rank them? HARD!!!; 5) allow missing data?; 7) allow model selection on traits.coefs using SSVS, with a possible extension of hypparams 
##############
#   rm(list = ls())
#   library(R2jags); 
#   library(mvtnorm); 
#   library(mvabund); 
#   library(coda); 
#   library(MASS)
#   library(fishMod)
#   source("auxilaryfunctions.R")

# family = "negative.binomial"; num.lv = 2; row.eff = "none"; n.burnin = 10000; n.iteration = 40000; n.thin = 30; save.model = TRUE; calc.ics = TRUE; trial.size <- 1; seed <- 123; hypparams = c(50,50,20,50); ssvs.index <- -1; do.fit = TRUE; model.name = NULL

#library(ade4)
#data(dunedata)
#y = dunedata$veg+1 ## Shift levels up to start at 1
#X = model.matrix(~A1 + factor(use) - 1, data = dunedata$envir)
#family = rep("ordinal",30)
#num.lv = 2; row.eff = "none"; n.burnin = 10000; n.iteration = 40000; n.thin = 30; save.model = TRUE; seed = 7; calc.ics = TRUE; trial.size = NULL; num.ord.levels <- 5; hypparams = c(50,20,50,50); 
#X <- matrix(rnorm(30*4),30,4)
#true.beta <- cbind(matrix(rnorm(length(family)*(ncol(X)+1)),length(family),ncol(X)+1),NA); 
#true.beta[nrow(true.beta),1] <- -sum(true.beta[-nrow(true.beta),1])
#true.ordinal <- seq(-0.5,0.5,length=num.ord.levels-1)
#y <- create.life(lv.coefs = true.beta[,c(1,ncol(true.beta))], X = X, X.coefs = true.beta[,-c(1,ncol(true.beta))], family = family, cutoffs = true.ordinal)

# n = 60; p <- 30
# X <- matrix(rnorm(n*2),n,2); beta <- cbind(matrix(rnorm(p*3),p,3),runif(p,0,5)); true.power <- 1.6
# mu <- exp(cbind(1,X)%*%t(beta[,1:3]))
# y <- matrix(NA,n,p)
# for(j in 1:ncol(y)) { y[,j] <- rTweedie(nrow(y), mu = mu[,j], phi = beta[j,4], p = true.power) }
# family = "tweedie"
# num.lv = 0; row.eff = "none"; n.burnin = 10000; n.iteration = 40000; n.thin = 30; save.model = TRUE; seed = 1; calc.ics = TRUE; trial.size = NULL; hypparams = c(50,20,50,50); 
# 
# library(FD)
# data(tussock)
# y <- tussock$trait[,c("height","LDMC","leafN","leafsize","SLA","seedmass","clonality","resprouting","lifespan")]
# y$LDMC <- y$LDMC/1000 ## change to g/g
# y$leafsize <- y$leafsize/100 ## change to cm^2
# y[,"resprouting"] <- as.numeric(y[,"resprouting"])-1 ## 0 = no; 1 = yes
# levels(y[,"clonality"]) = c(3,2,1)
# y[,"clonality"] <- as.numeric(levels(y[,"clonality"]))[y[,"clonality"]]
# levels(y[,"lifespan"]) = c(0,0,1)
# y[,"lifespan"] <- as.numeric(levels(y[,"lifespan"]))[y[,"lifespan"]]
# y <- y[-which(is.na(rowSums(y))),]
# 
# family = c("lnormal","normal","normal","lnormal","lnormal","lnormal","ordinal","binomial","binomial")
# num.lv = 2; row.eff = "none"; n.burnin = 10000; n.iteration = 40000; n.thin = 30; save.model = TRUE; seed = 123; calc.ics = TRUE; trial.size = 1; hypparams = c(100,100,100,100); X <- NULL

# n = 30; s <- 30; num.multinom.levels <- 4
# X <- matrix(rnorm(n*2),n,2)
# X.coefs <- rbind(matrix(0,s-1,2),c(1,2))
# X.multinom.coefs <- array(NA,dim=c(s-1,2,num.multinom.levels))
# for(k in 1:num.multinom.levels) { X.multinom.coefs[,,k] <- rnorm((s-1)*2) }
# row.coefs <- runif(n)
# lv.coefs <- cbind(matrix(runif(s,-3,-1),s,1),2)
# family = c(rep("multinom",s-1),"normal")
# num.lv = 2; row.eff = "fixed"; n.burnin = 10000; n.iteration = 40000; n.thin = 30; save.model = TRUE; seed = 1; calc.ics = TRUE; trial.size = 1; hypparams = c(50,20,50,50); 
# y <- create.life(lv.coefs = lv.coefs, X = X, X.coefs = X.coefs, X.multinom.coefs = X.multinom.coefs, family = family, row.coefs = row.coefs)

boral <- function(y, ...) UseMethod("boral")

## Model is g(mu_{ij}) = row + beta0 + LV_i*theta_j + X_i*beta_j
boral.default <- function (y, X = NULL, traits = NULL, which.traits = NULL, family, trial.size = 1, num.lv = 0, row.eff = "none", n.burnin = 10000, n.iteration = 40000, n.thin = 30, save.model = FALSE, seed = 123, calc.ics = TRUE, hypparams = c(100, 20, 100, 50), ssvs.index = -1, do.fit = TRUE, model.name = NULL, ...) {

	if(is.null(dim(y))) { 
		cat("Converting y into a one column matrix.\n"); y <- matrix(y, ncol = 1) }
	if(!is.null(X) & is.null(dim(X))) { 
		cat("Converting X into a one column matrix\n"); X <- matrix(X, ncol = 1) }
	if(!is.null(traits) & is.null(dim(traits))) { 
		cat("Converting traits into a one column matrix\n"); traits <- matrix(traits, ncol = 1) }
	if(length(hypparams) != 4) { 
		stop("hypparams must be a vector of four elements. Please see boral help file as to what the elements correspond to.\n") }
 	if(!is.null(which.traits)) { 
		print("Current version of boral ignores does not support the ssvs.index argument when which.traits is supplied. Sorry!"); ssvs.index <- -1 }
	if(!is.null(X)) { if(!is.matrix(X)) X <- as.matrix(X) }
	if(!is.null(X)) { if(any(apply(X,2,function(x) all(x == 1)))) { 
		stop("No intercept column should be included in X") } }
	if(!is.null(traits)) { if(!is.matrix(traits)) traits <- as.matrix(traits) }
    
    
	if(num.lv == 1) warning("We won't stop you, but one latent variable is unlikely to be successful in capturing between column correlation!")
	if(num.lv > 5) warning("We won't stop you, but please consider if you really want more than five latent variables in the model!")
	
	
	if(length(family) != ncol(y) & length(family) != 1) { stop("Number of elements in family must either one or the # of columns in y") }
	if(length(family) == 1) family <- rep(family, ncol(y))
	if(!all(family %in% c("negative.binomial", "poisson", "binomial", "normal", "lnormal", "tweedie", "ordinal", "exponential", "gamma", "beta"))) stop("At least one of the elements in family is not compatible with current version of boral...sorry!")
	if(any(family == "ordinal")) {
		if(sum(y[, family == "ordinal"] == 0) > 0) stop("For ordinal data, please shift minimum level to 1.")
		print("It is assumed all ordinal columns have the same number of levels -- please see help file as to the motivation behind this.")
		print("boral may take a ridiculously long time to fit ordinal data models. Apolgoies in advance!") 
		if(!is.null(traits)) stop("Current version of boral does not allow traits for ordinal responses. Sorry!")
		}
	
	
	if(row.eff == FALSE) row.eff <- "none"; 
	if(row.eff == TRUE) row.eff <- "fixed"
	if(!(row.eff %in% c("none", "fixed", "random"))) stop("row.eff must be one of none/fixed/random.")
	
	if(!is.null(X)) { num.X <- ncol(X) } else { num.X <- 0 }
	if(!is.null(traits)) { num.traits <- ncol(traits) } else { num.traits <- 0 }
	
	if(num.X == 0 & num.traits > 0) 
		stop("num.traits > 0 suggests traits are to be regressed against covariates X, so please supply X.") 
 	if(num.traits > 0 & is.null(which.traits)) 
		stop("If traits are supplied, then please also supply which.traits to inform what traits are regressed against which covariates.") 
 	if(!is.null(which.traits) & ((num.X+1) != length(which.traits))) 
		stop("which.traits should have equal to 1+length(ncol(X))") 
 	if(!is.null(which.traits)) { if(any(sapply(which.traits,length) > num.traits)) stop("Each element in the list which.traits should have at most ncol(traits) elements.") }
 	#if(is.null(which.traits)) { which.traits <- vector("list",num.X+1); for(k in 1:length(num.X+1)) which.traits[[k]] <- 0 } 

	
	if(!(length(ssvs.index) %in% c(1, ncol(X)))) 
		stop("Number of elements in ssvs.index must either be one or the # of columns in X.")
	if(length(ssvs.index) == 1 & num.X > 0) ssvs.index <- rep(ssvs.index, ncol(X))
	if(any(ssvs.index < -1)) 
		stop("Elements of ssvs.index can only take values in -1, 0, or any positive integer; please see help file for guide.")
	
	
	if(any(family == "binomial") & !(length(trial.size) %in% c(1, length(family)))) 
        stop("trial.size needs to be specified if any columns are binomially distributed; can either be a single element or a vector equal to the # of columns in y. The latter will assume the specified trial size for all rows labelled binomial in the family argument.")
	if(any(family == "binomial") & length(trial.size) == 1) {
		complete.trial.size <- rep(0, ncol(y))
		complete.trial.size[which(family == "binomial")] <- trial.size }
	if(any(family == "binomial") & length(trial.size) == ncol(y)) { complete.trial.size <- trial.size }
	if(all(family != "binomial")) { complete.trial.size <- rep(0, ncol(y)) }
	if(all(family == "binomial") & all(complete.trial.size == 1)) { family <- rep("bernoulli",ncol(y)) }
	
	
	if(all(family != "ordinal")) { num.ord.levels <- 0; }
	if(any(family == "ordinal")) { num.ord.levels <- max(y[, family == "ordinal"]); }
	if(all(family != "multinom")) { num.multinom.levels <- 0; index.multinom.cols <- NULL }
# 	if(any(family == "multinom")) { 
# 		num.multinom.levels <- apply(y[, family == "multinom"], 2, max)
# 		index.multinom.cols <- which(family == "multinom") 
# 		}

		
	n <- nrow(y); p <- ncol(y)
 	#n.chains <- 1; ## Run one chain only to avoid arbitrary rotation problems

 	
 	if(num.lv > 0) 
		make.jagsboralmodel(family, num.X, num.traits, which.traits, row.eff, complete.trial.size, n, p, hypparams, ssvs.index, model.name)
	if(num.lv == 0)  
		make.jagsboralnullmodel(family, num.X, num.traits, which.traits, row.eff, complete.trial.size, n, p, hypparams, ssvs.index, model.name)
 	if(!do.fit) { 
		cat("JAGS model file created only. Thank you, come again!\n")
		return() }

		
	jags.data <- list("y", "n", "p", "num.lv", "num.X", "num.traits", "num.ord.levels", "num.multinom.levels")
	if(num.X > 0) jags.data <- c("X", jags.data)
	if(num.traits > 0) jags.data <- c("traits", jags.data)
	if(any(family == "ordinal")) { ones <- matrix(1, n, p); jags.data <- c(jags.data, "ones") }
	
	
	jags.params <- c("all.params")
	if(num.lv > 0) jags.params <- c(jags.params, "lvs")
	if(row.eff != "none") jags.params <- c(jags.params, "row.params")
	if(row.eff == "random") jags.params <- c(jags.params, "row.ranef.sigma")
	if(num.X > 0 & any(family != "multinom")) jags.params <- c(jags.params, "X.params")
	#if(num.X > 0 & any(family == "multinom")) jags.params <- c(jags.params, "X.multinom.params")
	if(num.traits > 0) jags.params <- c(jags.params, "traits.params", "sigma.trait")
	if(any(family == "tweedie")) jags.params <- c(jags.params, "powerparam")
	if(any(family == "ordinal")) jags.params <- c(jags.params, "alpha")
	if(any(ssvs.index == 0)) jags.params <- c(jags.params, paste("probindX", which(ssvs.index == 0), sep = ""))
	if(any(ssvs.index > 0)) jags.params <- c(jags.params, paste("probGpX", unique(ssvs.index[ssvs.index > 0]), sep = ""))
	
	
	jags.inits <- function() {
		initial.list <- list()
		#if(num.lv > 0) initial.list$lvs <- matrix(rnorm(n*num.lv,0,0.1),n,num.lv)
#  		if(num.lv > 0) { 
#  			initial.list$all.params <- matrix(0,p,num.lv+1) 
# 			if(!all(family %in% c("poisson", "binomial", "ordinal", "multinom", "exponential","bernoulli"))) { initial.list$all.params <- cbind(initial.list$all.params,0.01) }
# 			}
		if(any(family %in% "tweedie")) initial.list$numfish = matrix(1, n, sum(family=="tweedie"))
		if(any(family %in% "ordinal")) initial.list$alpha0 <- seq(-1, 1, length = num.ord.levels - 1)
		
		if(all(family %in% "bernoulli")) {
			Tau <- rWishart(1,p+1,diag(p))[,,1]
			Sigma <- solve(Tau)
			Z <- abs(t(rmvnorm(n,rep(0,p),Sigma)))
			Z <- ifelse(as.matrix(y), Z, -1 * Z)
			initial.list$Z <- Z }
			
		return(initial.list)
		}

	set.seed(seed)
	actual.filename <- model.name
	if(is.null(actual.filename)) actual.filename <- "jagsboralmodel.txt"


	## The fit! ##
	n.chains <- 1
	jagsfit <- try(suppressWarnings(jags(data = jags.data, inits = jags.inits, jags.params, model.file = actual.filename, n.iter = n.iteration, n.burnin = n.burnin, n.chains = n.chains, n.thin = n.thin)),silent=TRUE)
	#jagsfit <- jags(data = jags.data, inits = jags.inits, parameters = jags.params, model.file = actual.filename)
    
	#print(jagsfit)
	if(inherits(jagsfit,"try-error")) {
		lookfornberror <- grep("Slicer stuck at value with infinite density", jagsfit[[1]])
		if(any(family == "negative.binomial") & length(lookfornberror) == 1) { 
			cat("MCMC fitting through JAGS failed. This is likely due to the prior on the dispersion (size) parameter of the negative binomial distribution been too uninformative. Please consider a tougher prior or switch to a Poisson family for those response that don't appear to actually be overdispersed. The error below informs you explicity which column of y the MCMC sampling ran into trouble:\n\n")
			print(jagsfit) }

		else {
			cat("MCMC fitting through JAGS failed:\n")
			print(jagsfit) }

		cat("boral fit failed...Exiting. Sorry!\n") 
		return()
		}
    
    
	## Format into big matrix; 
	fit.mcmcBase <- jagsfit$BUGSoutput
	fit.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = n.thin) ## Thanks to Guilliaume Blanchet for this!
	if(n.chains == 1) {
		combined.fit.mcmc <- fit.mcmc
# 		gdiagnostic <- geweke.diag(combined.fit.mcmc)
		}
# 	if(n.chains > 1) {
# 		combined.fit.mcmc <- as.mcmc(fit.mcmcBase$sims.matrix)
# 
# 		fit.rhats <- (rhats(jagsfit, asc = FALSE))
# 		make.rhatslist <- list(lv.coefs = matrix(fit.rhats[grep("all.params", rownames(fit.rhats)),], nrow = p))
# 		fit.rhats <- fit.rhats[-grep("all.params",rownames(fit.rhats)),]
# 		if(num.lv > 0) { make.rhatslist$lvs <- fit.rhats[grep("lvs",names(fit.rhats))] } 
# 		#print(fit.rhats)
# 		rownames(make.rhatslist$lv.coefs) <- colnames(y); colnames(make.rhatslist$lv.coefs) <- NULL
# 	
# 		if(row.eff != "none") {
# 			make.rhatslist$row.coefs <- fit.rhats[grep("row.params", names(fit.rhats))]		
# 			names(make.rhatslist$row.coefs) <- rownames(y) }
# 		if(row.eff == "random") {
# 			make.rhatslist$row.ranef <- fit.rhats[grep("row.ranef.sigma", names(fit.rhats))]
# 			names(make.rhatslist$row.ranef) <- c("Row random effects sigma") }
# 		
# 		if(num.X > 0) {
# 			make.rhatslist$X.coefs <- matrix(fit.rhats[grep("X.params", names(fit.rhats))], nrow = p)
# 			rownames(make.rhatslist$X.coefs) <- colnames(y); colnames(make.rhatslist$X.coefs) <- colnames(X) }
# 		
# 		if(any(family == "ordinal")) {
# 			make.rhatslist$cutoffs <- fit.rhats[grep("alpha", names(fit.rhats))]
# 			names(make.rhatslist$cutoffs) <- paste(1:(num.ord.levels - 1), "|", 2:num.ord.levels, sep = "") }
# 
# 		if(any(family == "tweedie")) {
# 			make.rhatslist$powerparam <- fit.rhats[grep("powerparam", names(fit.rhats))]
# 			names(make.rhatslist$powerparam) <- "Common power parameter" }
# 	
# 		#print(make.rhatslist)
# 		exceed.rhatcutoff <- sum(sapply(make.rhatslist, function(x) sum(x > rhat.cutoff)))
# 		cat("There were", exceed.rhatcutoff, "(", 100*exceed.rhatcutoff/sum(sapply(make.rhatslist,length)), "%) parameters whose Rhat exceeded the prespecified cutoff of", rhat.cutoff, "\n")		
# 		rm(fit.rhats)
# 		}
	rm(fit.mcmc)
    
    
#   	## Flip dispersion parameters?
#  	sel.thetas <- grep("all.params", colnames(combined.fit.mcmc))
#  	sel.thetas2 <- as.numeric(sel.thetas[(length(sel.thetas) - p + 1):length(sel.thetas)])		
#  	if(any(family %in% c("negative.binomial"))) {
#  		combined.fit.mcmc[, sel.thetas2[family %in% c("negative.binomial")]] <- 1/combined.fit.mcmc[, sel.thetas2[family %in% c("negative.binomial")]] 
#  		}

		
#   	## For any multinomial columns, set the corresponding rows in X.coefs to zero
# 	if(any(family == "multinom") & num.X > 0) {
# 		for(k in index.multinom.cols) {
# 			sel.multinom.col <- grep(paste("X.params\\[", k, ",+", sep = ""), colnames(combined.fit.mcmc))
# 			combined.fit.mcmc[, sel.multinom.col] <- 0 }
# 		}

		
 	## Make output beautiful
	if(is.null(colnames(y))) colnames(y) <- 1:ncol(y); if(is.null(rownames(y))) rownames(y) <- 1:nrow(y)
	if(num.X > 0) { if(is.null(colnames(X))) colnames(X) <- 1:ncol(X); if(is.null(rownames(X))) rownames(X) <- 1:nrow(X) }
	if(num.traits > 0) { 
		if(is.null(colnames(traits))) colnames(X) <- 1:ncol(traits); if(is.null(rownames(traits))) rownames(X) <- 1:nrow(traits) }

		
	out.fit <- list(lv.coefs.median = matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, median), nrow = p), lv.coefs.iqr = matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, IQR), nrow = p), lv.coefs.mean = matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, mean), nrow = p), lv.coefs.sd = matrix(apply(combined.fit.mcmc[, grep("all.params", colnames(combined.fit.mcmc))], 2, sd), nrow = p))
	rownames(out.fit$lv.coefs.median) <- rownames(out.fit$lv.coefs.iqr) <- rownames(out.fit$lv.coefs.mean) <- rownames(out.fit$lv.coefs.sd) <- colnames(y)

	
	if(num.lv > 0) { 
		out.fit$lv.median = matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, median), nrow = n)
		out.fit$lv.iqr = matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, IQR), nrow = n)
		out.fit$lv.mean = matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, mean), nrow = n)
		out.fit$lv.sd = matrix(apply(combined.fit.mcmc[, grep("lvs", colnames(combined.fit.mcmc))], 2, sd), nrow = n)
		rownames(out.fit$lv.median) <- rownames(out.fit$lv.iqr) <- rownames(out.fit$lv.mean) <- rownames(out.fit$lv.sd) <- rownames(y)
		colnames(out.fit$lv.median) <- colnames(out.fit$lv.iqr) <- colnames(out.fit$lv.mean) <- colnames(out.fit$lv.sd) <- paste("LV", 1:num.lv, sep = "")
		
		if(ncol(out.fit$lv.coefs.median) == (num.lv+2)) colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.sd) <- c("beta0", paste("theta", 1:num.lv, sep = ""), "Dispersion") 
		if(ncol(out.fit$lv.coefs.median) == (num.lv+1)) colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.sd) <- c("beta0", paste("theta", 1:num.lv, sep = "")) 
		}
	if(num.lv == 0) {
		if(ncol(out.fit$lv.coefs.median) == 2) colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.sd) <- c("beta0", "Dispersion") 
		if(ncol(out.fit$lv.coefs.median) == 1) colnames(out.fit$lv.coefs.median) <- colnames(out.fit$lv.coefs.iqr) <- colnames(out.fit$lv.coefs.mean) <- colnames(out.fit$lv.coefs.sd) <- c("beta0") 
		}	
	
	if(row.eff != "none") {
		out.fit$row.coefs.median <- apply(combined.fit.mcmc[, grep("row.params", colnames(combined.fit.mcmc))], 2, median)
		out.fit$row.coefs.iqr <- apply(combined.fit.mcmc[, grep("row.params", colnames(combined.fit.mcmc))], 2, IQR)
		out.fit$row.coefs.mean <- apply(combined.fit.mcmc[, grep("row.params", colnames(combined.fit.mcmc))], 2, mean)
		out.fit$row.coefs.sd <- apply(combined.fit.mcmc[, grep("row.params", colnames(combined.fit.mcmc))], 2, sd)
		
		names(out.fit$row.coefs.median) <- names(out.fit$row.coefs.iqr) <- names(out.fit$row.coefs.mean) <- names(out.fit$row.coefs.sd) <- rownames(y)
	
		if(row.eff == "random") {
			out.fit$row.sigma.median <- median(combined.fit.mcmc[, grep("row.ranef.sigma", colnames(combined.fit.mcmc))])
			out.fit$row.sigma.iqr <- IQR(combined.fit.mcmc[, grep("row.ranef.sigma", colnames(combined.fit.mcmc))])
			out.fit$row.sigma.mean <- mean(combined.fit.mcmc[, grep("row.ranef.sigma", colnames(combined.fit.mcmc))])
			out.fit$row.sigma.sd <- sd(combined.fit.mcmc[, grep("row.ranef.sigma", colnames(combined.fit.mcmc))])
            
			names(out.fit$row.sigma.median) <- names(out.fit$row.sigma.iqr) <- names(out.fit$row.sigma.mean) <- names(out.fit$row.sigma.sd) <- c("Row random effects sigma") }
		}

		
	if(num.X > 0) {
		out.fit$X.coefs.median <- matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, median), nrow = p)
		out.fit$X.coefs.iqr <- matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, IQR), nrow = p)
		out.fit$X.coefs.mean <- matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, mean), nrow = p)
		out.fit$X.coefs.sd <- matrix(apply(combined.fit.mcmc[, grep("X.params", colnames(combined.fit.mcmc))], 2, sd), nrow = p)
		rownames(out.fit$X.coefs.median) <- rownames(out.fit$X.coefs.iqr) <- rownames(out.fit$X.coefs.mean) <- rownames(out.fit$X.coefs.sd) <- colnames(y)
		colnames(out.fit$X.coefs.median) <- colnames(out.fit$X.coefs.iqr) <- colnames(out.fit$X.coefs.mean) <- colnames(out.fit$X.coefs.sd) <- colnames(X)
		
		if(any(ssvs.index == 0)) {
			out.fit$ssvs.indcoefs.mean <- matrix(apply(combined.fit.mcmc[, grep("probindX", colnames(combined.fit.mcmc))], 
                2, mean), nrow = p)
			rownames(out.fit$ssvs.indcoefs.mean) <- colnames(y)
			colnames(out.fit$ssvs.indcoefs.mean) <- colnames(X)[which(ssvs.index == 0)]
			out.fit$ssvs.indcoefs.sd <- matrix(apply(combined.fit.mcmc[, grep("probindX", colnames(combined.fit.mcmc))], 2, sd), nrow = p)
			rownames(out.fit$ssvs.indcoefs.sd) <- colnames(y)
			colnames(out.fit$ssvs.indcoefs.sd) <- colnames(X)[which(ssvs.index == 0)]
			}
		if(any(ssvs.index > 0)) {
			out.fit$ssvs.gpcoefs.mean <- apply(as.matrix(combined.fit.mcmc[, grep("probGpX", colnames(combined.fit.mcmc))]), 2, mean)
			names(out.fit$ssvs.gpcoefs.mean) <- paste("Gp", unique(ssvs.index[ssvs.index > 0]), sep = "")
			out.fit$ssvs.gpcoefs.sd <- apply(as.matrix(combined.fit.mcmc[, grep("probGpX", colnames(combined.fit.mcmc))]), 2, sd)
			names(out.fit$ssvs.gpcoefs.sd) <- paste("Gp", unique(ssvs.index[ssvs.index > 0]), sep = "")
			}
		}
		
		
	if(num.traits > 0) {
		out.fit$traits.coefs.median <- cbind(matrix(apply(combined.fit.mcmc[, grep("traits.params", colnames(combined.fit.mcmc))], 2, median), nrow = num.X+1),apply(combined.fit.mcmc[, grep("sigma.trait", colnames(combined.fit.mcmc))], 2, median))
		out.fit$traits.coefs.iqr <- cbind(matrix(apply(combined.fit.mcmc[, grep("traits.params", colnames(combined.fit.mcmc))], 2, IQR), nrow = num.X+1),apply(combined.fit.mcmc[, grep("sigma.trait", colnames(combined.fit.mcmc))], 2, IQR))
		out.fit$traits.coefs.mean <- cbind(matrix(apply(combined.fit.mcmc[, grep("traits.params", colnames(combined.fit.mcmc))], 2, mean), nrow = num.X+1),apply(combined.fit.mcmc[, grep("sigma.trait", colnames(combined.fit.mcmc))], 2, mean))
		out.fit$traits.coefs.sd <- cbind(matrix(apply(combined.fit.mcmc[, grep("traits.params", colnames(combined.fit.mcmc))], 2, sd), nrow = num.X+1),apply(combined.fit.mcmc[, grep("sigma.trait", colnames(combined.fit.mcmc))], 2, sd))

		rownames(out.fit$traits.coefs.median) <- rownames(out.fit$traits.coefs.iqr) <- rownames(out.fit$traits.coefs.mean) <- rownames(out.fit$traits.coefs.sd) <- c("beta0",colnames(X))
		colnames(out.fit$traits.coefs.median) <- colnames(out.fit$traits.coefs.iqr) <- colnames(out.fit$traits.coefs.mean) <- colnames(out.fit$traits.coefs.sd) <- c(colnames(traits),"sigma")
		}

#   	if(num.X > 0 & any(family == "multinom")) {
#   		out.fit$X.multinom.coefs.median <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,median),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
#   		out.fit$X.multinom.coefs.iqr <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,IQR),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
#   		out.fit$X.multinom.coefs.mean <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,mean),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
#   		out.fit$X.multinom.coefs.sd <- array(apply(combined.fit.mcmc[,grep("X.multinom.params", colnames(combined.fit.mcmc))],2,sd),dim=c(length(index.multinom.cols),num.X,num.multinom.levels))
#   
#   		dimnames(out.fit$X.multinom.coefs.median) <- dimnames(out.fit$X.multinom.coefs.iqr) <- dimnames(out.fit$X.multinom.coefs.mean) <- dimnames(out.fit$X.multinom.coefs.sd) <- list("1" = index.multinom.cols, "2" = colnames(X), "level" = 1:num.multinom.levels)
#   		}

	if(any(family == "ordinal")) {
		out.fit$cutoffs.median <- apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, median)
		out.fit$cutoffs.iqr <- apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, IQR)
		out.fit$cutoffs.mean <- apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, mean)
		out.fit$cutoffs.sd <- apply(combined.fit.mcmc[, grep("alpha", colnames(combined.fit.mcmc))], 2, sd)
		names(out.fit$cutoffs.median) <- names(out.fit$cutoffs.iqr) <- names(out.fit$cutoffs.mean) <- names(out.fit$cutoffs.sd) <- paste(1:(num.ord.levels - 1), "|", 2:num.ord.levels, sep = "") }

	if(any(family == "tweedie")) {
		out.fit$powerparam.median <- median(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))])
		out.fit$powerparam.iqr <- IQR(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))])
		out.fit$powerparam.mean <- mean(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))])
		out.fit$powerparam.sd <- sd(combined.fit.mcmc[, grep("powerparam", colnames(combined.fit.mcmc))])
		names(out.fit$powerparam.median) <- names(out.fit$powerparam.iqr) <- names(out.fit$powerparam.mean) <- names(out.fit$powerparam.sd) <- "Common power parameter"
		}

	#print(out.fit$lv.coefs.mean)
	get.hpds <- get.hpdintervals(y, X, traits, combined.fit.mcmc, num.lv)
	out.fit$hpdintervals <- get.hpds
	if(calc.ics) {
		cat("Calculating Information criteria\n")
		if(num.traits > 0) 
			print("Please note that AIC and BIC at post median are currently calculated without taking into account the traits, i.e. that the column-specific coefficients are now random effects!")
		get.ics <- get.measures(y, X, family, complete.trial.size, row.eff, num.lv, combined.fit.mcmc, more.measures = FALSE)
		ics <- c(get.dic(jagsfit), get.ics$waic, get.ics$eaic, get.ics$ebic, get.ics$aic.median, get.ics$bic.median)#, get.ics$comp.lm)
		names(ics) <- c("Conditional DIC", "WAIC", "EAIC", "EBIC", "AIC at post. median", "BIC at post. median")#, "Compound L-M at post. median")
		out.fit$ics <- ics
		}
			
	if(save.model) { out.fit$jags.model <- jagsfit }

	out.fit$call <- match.call()
	out.fit$n <- n; out.fit$p <- p
	out.fit$X <- X
	out.fit$traits <- traits
	out.fit$y <- y

	out.fit$family <- family; if(all(family == "bernoulli")) out.fit$family <- rep("binomial",p)
	out.fit$num.lv <- num.lv
	out.fit$num.X <- num.X; out.fit$num.traits <- num.traits
	out.fit$which.traits <- which.traits
	out.fit$row.eff <- row.eff
	out.fit$calc.ics <- calc.ics
	out.fit$trial.size <- complete.trial.size
	out.fit$hypparams <- hypparams
	out.fit$ssvs.index <- ssvs.index
	out.fit$num.ord.levels <- num.ord.levels
	out.fit$n.burnin <- n.burnin; out.fit$n.thin <- n.thin; out.fit$n.iteration <- n.iteration; #out.fit$n.chain <- out.fit$n.chains; 
	#if(n.chains == 1) out.fit$geweke.diag <- gdiagnostic
	
	class(out.fit) <- "boral"
	if(!save.model) { if(file.exists(actual.filename)) file.remove(actual.filename) }

	return(out.fit) }
 	

 	
################	
lvsplot <- function(x, jitter = FALSE, a = 1, biplot = TRUE, ind.spp = NULL, alpha = 0.5, main = NULL, est = "median", ...) {
 	if(x$num.lv > 2) stop("Manual plotting required for plotting beyond 2 latent variables.")
 	if(x$num.lv == 0) stop("No latent variables to plot.")
 
 	n <- nrow(x$lv.median); p <- nrow(x$lv.coefs.median)
 	if(!is.null(ind.spp)) { if(ind.spp > p) { ind.spp <- p } }
	if(biplot == TRUE & !is.null(ind.spp)) { 
		cat("Only the first", ind.spp, "`most important' latent variable coefficients included in biplot\n") }
	if(biplot == TRUE & is.null(ind.spp)) { 
		ind.spp <- p; cat("All latent variable coefficients included in biplot\n") }

 	par(cex = a, cex.axis = a, cex.lab = a+0.5, mar = c(5,5,3,1), mfrow = c(1,1), las = 1, cex.main = a+0.5, ...) 
 
 	if(x$num.lv == 1) {
		choose.x <- x$lv.median; 
		main <- "Plot of latent variable posterior medians"
		if(est == "mean") { choose.x <- x$lv.mean; main <- "Plot of latent variable posterior means" }
		plot(1:n, choose.x, xlab = "Row index", ylab = "Latent variable 1", main = main, cex = 1.2*a, type = "n", ...)
		text(x = 1:n, y = x$lv.median, label = 1:n, cex = 1.2*a)
		}


 	if(x$num.lv == 2) {
# 		## Scale by L2norms
# 		x$lv.median2 <- scale(x$lv.median,center=TRUE,scale=FALSE) #*matrix(1/sqrt(colSums(x$lv.median^2)),n,2,byrow=TRUE)
# 		x$lv.coefs.median2 <- scale(x$lv.coefs.median[,2:3]*matrix(sqrt(colSums(x$lv.median2^2))/sqrt(colSums(x$lv.coefs.median[,2:3]^2)),p,2,byrow=TRUE),center=TRUE,scale=FALSE) 

#  		get.lv.norms <- sqrt(colSums(x$lv.median^2))
#  		get.lv.coefs.norms <- sqrt(colSums(x$lv.coefs.median[,2:3]^2))
# 		D <- diag(x=get.lv.norms,2,2)*diag(x=get.lv.coefs.norms,2,2)
#  		x$lv.median2 <- (x$lv.median*matrix(1/get.lv.norms,x$n,2,byrow=TRUE))%*%D^alpha
#  		x$lv.coefs.median2 <- (x$lv.coefs.median[,2:3]*matrix(1/get.lv.coefs.norms,x$p,2,byrow=TRUE))%*%D^(1-alpha)
 		
   		testcov <- x$lv.median%*%t(x$lv.coefs.median[,2:3])
		if(est == "mean") { testcov <- x$lv.mean%*%t(x$lv.coefs.mean[,2:3]) }

		do.svd <- svd(testcov,x$num.lv,x$num.lv)   		
   		choose.lvs <- scale(do.svd$u*matrix(do.svd$d[1:x$num.lv]^alpha,nrow=x$n,ncol=2,byrow=T),center=T, scale = F)
   		choose.lv.coefs <- scale(do.svd$v*matrix(do.svd$d[1:x$num.lv]^(1-alpha),nrow=x$p,ncol=2,byrow=T),center=T, scale = F)
   		
#  		## Scale by L2norms
# 		x$lv.mean2 <- scale(x$lv.mean,center=TRUE,scale=FALSE) #*matrix(sqrt(colSums(x$lv.coefs.mean[,2:3]^2))/sqrt(colSums(x$lv.mean^2)),n,2,byrow=TRUE) 
# 		x$lv.coefs.mean2 <- scale(x$lv.coefs.mean[,2:3]*matrix(sqrt(colSums(x$lv.mean2^2))/sqrt(colSums(x$lv.coefs.mean[,2:3]^2)),p,2,byrow=TRUE),center=FALSE,scale=FALSE) 

		if(!biplot) {
			if(is.null(main) & est == "median") { main = "Plot of latent variable posterior medians" }
			if(is.null(main) & est == "mean") { main = "Plot of latent variable posterior means" }
			plot(choose.lvs, xlab = "Latent variable 1", ylab = "Latent variable 2", main = main, cex = 1.2*a, type = "n", ...)
			if(!jitter) text(choose.lvs, label = 1:n, cex = 1.2*a)
			if(jitter) text(jitter(choose.lvs[,1]), jitter(choose.lvs[,2]), label = 1:n, cex = 1.2*a)
			}

		if(biplot) {
			if(is.null(main) & est == "median") { main = "Biplot of latent variable posterior medians" }
			if(is.null(main) & est == "mean") { main = "Biplot of latent variable posterior means" }
			largest.lnorms <- order(rowSums(choose.lv.coefs^2),decreasing=TRUE)[1:ind.spp]	
			
			plot(rbind(choose.lvs,choose.lv.coefs), xlab = "Latent variable 1", ylab = "Latent variable 2", main = main, cex = a, type = "n", xlim = 1.1*range(rbind(choose.lvs,choose.lv.coefs)[,1]), ylim = 1.1*range(rbind(choose.lvs,choose.lv.coefs)[,2]))
			if(!jitter) text(choose.lvs, label = 1:n, cex = 1.2*a)
			if(jitter) text(jitter(choose.lvs[,1]), jitter(choose.lvs[,2]), label = 1:n, cex = 1.2*a)
			text(choose.lv.coefs[largest.lnorms,], label = rownames(x$lv.coefs.mean[largest.lnorms,]), col = "red", cex = 0.9*a)	
			}
 		}	

 	par(mfrow = c(1,1), las = 0, cex = 1, cex.axis = 1, cex.lab = 1, ask = FALSE, cex.main = 1.2, mar = c(5.1,4.1,4.1,2.1)) 	
 	}


print.boral <- function(x, ...) {
 	cat("Call:\n")
 	print(x$call)
 	cat("\n")
 	cat("Response matrix attributes\n \t# of rows:", x$n, "\n\t# of columns:", x$p, "\n") 
 	cat("Model attributes\n \tColumn families:", x$family, "\n\t# of latent variables:", x$num.lv, "\n\tRow effect included (none/fixed/random)?", x$row.eff, "\n") 
 	if(any(x$family == "binomial")) cat("Trial sizes used (columns with binomial families):", x$trial.size,"\n")
 	if(any(x$family == "ordinal")) cat("Number of levels for ordinal data:", x$num.ord.levels,"\n")
 	if(x$num.X > 0) cat("Model matrix with", x$num.X, "covariates also fitted\n\n")
 	if(x$num.traits > 0) cat("Trait matrix with", x$num.traits, "traits also included\n\n")
 	if(any(x$ssvs.index > -1)) cat("SSVS performed on covariates with indices", x$ssvs.index, "\n\n")
# 	cat("Output attributes\n")
# 	print(attributes(x))
 	}

 	
print.summary.boral <- function(x, ...) {
 	cat("Call:\n")
 	print(x$call)
 	cat("\n")
 	
 	if(x$est == "median") { cat("Median point estimates\n\n LV coefficients (thetas) and dispersion parameter (if applicable)\n"); print(x$coefficients); cat("\n") }
 	if(x$est == "mean") { cat("Mean point estimates\n\n LV coefficients (thetas) and dispersion parameter (if applicable)\n"); print(x$coefficients); cat("\n") }	
 	
 	if(!is.null(x$row.coefficients)) { cat("Row coefficients\n"); print(x$row.coefficients); cat("\n") }
 	if(!is.null(x$X.coefficients)) { cat("X coefficients (betas)\n"); print(x$X.coefficients); cat("\n") }
 	if(!is.null(x$X.multinom.coefficients)) { cat("There are also coefficients corresponding to multinomial columns which have not been printed"); }
 	if(!is.null(x$traits.coefficients)) { cat("Trait coefficients\n"); print(x$traits.coefficients); cat("\n") }
 	
 	if(any(x$family == "ordinal")) { cat("Proportional odds (Cumulative probit) regression intercepts\n"); print(x$cutoffs); cat("\n") }
 	if(any(x$family == "tweedie")) { cat("Tweedie power parameter \n"); print(x$powerparam); cat("\n") }
 
 	if(x$calc.ics) {
 		cat("DIC (pD = var(deviance)/2):", as.vector(unlist(x$ics[1])), "\n")
 		cat("WAIC:", as.vector(unlist(x$ics[2])), "\n")
 		cat("EAIC:", as.vector(unlist(x$ics[3])), "\n")
 		cat("EBIC:", as.vector(unlist(x$ics[4])), "\n")
 		cat("AIC at posterior median:", as.vector(unlist(x$ics[5])), "\n")
 		cat("BIC at posterior median:", as.vector(unlist(x$ics[6])), "\n")
 		#cat("Compound Laplace-Metropolis estimator at posterior median:", as.vector(unlist(x$ics[7])), "\n")
		}

	if(!is.null(x$ssvs.indcoefs.prob)) { cat("SSVS probabilities on individual coefficients\n"); print(x$ssvs.indcoefs.prob); cat("\n") }
	if(!is.null(x$ssvs.gpcoefs.prob)) { cat("SSVS probabilities on groups of coefficients\n"); print(x$ssvs.gpcoefs.prob); cat("\n") }			
	}	
		
		
summary.boral <- function(object, est = "median", ...) {
	if(est == "median") {
 		gather.output <- list(call = object$call, coefficients = round(object$lv.coefs.median,3))
 		if(object$row.eff != "none") gather.output$row.coefficients = round(object$row.coefs.median,3)
 		if(object$num.X > 0) gather.output$X.coefficients = round(object$X.coefs.median,3)
 		if(object$num.traits > 0) gather.output$traits.coefficients = round(object$traits.coefs.median,3)
 		if(any(object$family == "ordinal")) gather.output$cutoffs = round(object$cutoffs.median,3)
 		if(any(object$family == "tweedie")) gather.output$powerparam = round(object$powerparam.median,3)
 		if(!is.null(object$X.multinom.coefs.median)) gather.output$X.multinom.coefficients = round(object$X.multinom.coefs.median,3) }
 
 	if(est == "mean") {
 		gather.output <- list(call = object$call, coefficients = round(object$lv.coefs.mean,3))
 		if(object$row.eff != "none") gather.output$row.coefficients = round(object$row.coefs.mean,3)
 		if(object$num.X > 0) gather.output$X.coefficients = round(object$X.coefs.mean,3)
 		if(object$num.traits > 0) gather.output$traits.coefficients = round(object$traits.coefs.mean,3)
 		if(any(object$family == "ordinal")) gather.output$cutoffs = round(object$cutoffs.mean,3)
 		if(any(object$family == "tweedie")) gather.output$powerparam = round(object$powerparam.mean,3)
 		if(!is.null(object$X.multinom.coefs.mean)) gather.output$X.multinom.coefficients = round(object$X.multinom.coefs.mean,3) }
 
 
 	gather.output$est <- est
 	gather.output$calc.ics <- object$calc.ics
	gather.output$trial.size <- object$trial.size
 	gather.output$num.ord.levels <- object$num.ord.levels
 	gather.output$ssvs.index <- object$ssvs.index
 
 
	if(any(object$ssvs.index == 0)) gather.output$ssvs.indcoefs.prob <- round(object$ssvs.indcoefs.mean,3)
	if(any(object$ssvs.index > 0)) gather.output$ssvs.gpcoefs.prob <- round(object$ssvs.gpcoefs.mean,3) 

	if(object$calc.ics) gather.output$ics <- object$ics
 	class(gather.output) <- "summary.boral"
 	gather.output 
 	}
 			
 			
plot.boral <- function(x, est = "median", jitter = FALSE, a = 1, ...) {
 	if(all(x$family %in% c("ordinal","multinom"))) 
 		stop("Residuals are not defined, and therefore residual analysis cannot be performed, if all columns of y are ordinal.")
 	get.mus <- fitted.boral(x, est = est)$out
 	get.etas <- get.mus
 	get.ds.res <- ds.residuals(object = x, est = est)$residuals
	for(j in 1:ncol(x$y)) {
 		if(x$family[j] %in% c("binomial","beta")) get.etas[,j] <- log((get.mus[,j]+1e-5)/(1-get.mus[,j]+1e-5))
 		if(x$family[j] %in% c("poisson","lnormal","negative.binomial","tweedie","exponential")) get.etas[,j] <- log(get.mus[,j]+1e-5)
 		if(x$family[j] == "normal") get.etas[,j] <- (get.mus[,j]) }
 
	par(ask = TRUE, cex = a, mar = c(5,5,2,1), cex.lab = 0.8*a, cex.main = a, las = 1, ...) 
 	palette(rainbow(ncol(get.etas)))
 
 	matplot(get.etas, get.ds.res, ylab = "Dunn-Smyth Residuals", xlab = "Linear Predictors", type="n")
 	for(i in 1:ncol(get.etas)) { points(get.etas[,i], get.ds.res[,i], col=palette()[i], cex = a) }
 	abline(h=0, lty = 2, lwd = 2)
	
# 	matplot(get.mus, get.ds.res, ylab = "Dunn-Smyth Residuals", xlab = "Fitted Values", type="n")
# 	for(i in 1:ncol(get.mus)) { points(get.mus[,i], get.ds.res[,i], col=palette()[i]) }
# 	abline(h=0, lty = 2, lwd = 2)father

	matplot(get.ds.res, ylab = "Dunn-Smyth Residuals", xlab = "Row index",type="n", xaxt = "n")
 	axis(side = 1, at = 1:nrow(x$y), labels = rownames(x$lv.mean), cex.axis = 1)
 	for (i in 1:ncol(get.mus)) { points(seq(1,nrow(x$y)),get.ds.res[,i], col=palette()[i], cex = a) }
 	abline(0,0,lty=2)
 
 	matplot(t(get.ds.res), ylab = "Dunn-Smyth Residuals", xlab = "Column index", type="n", xaxt = "n")
 	axis(side = 1, at = 1:ncol(x$y), labels = rownames(x$coefs.mean), cex.axis = 1)
 	for(i in 1:ncol(get.mus)) { points(rep(i,nrow(get.etas)), get.ds.res[,i], col=palette()[i], cex = a) }
 	abline(h=0, lty = 2, lwd = 2)
 
	get.ds.res2 <- as.vector(unlist(get.ds.res))
 	qqnorm(get.ds.res2[is.finite(get.ds.res2)], main = "Normal Quantile Plot", cex.main = a*0.9, cex = a)
 	#qqline(get.ds.res2[is.finite(get.ds.res2)], lwd = a)
 	abline(0,1, lwd = a)
 	
 	palette("default")
 	par(mfrow = c(1,1), las = 0, cex = 1, cex.axis = 1, cex.lab = 1, ask = FALSE, cex.main = 1.2, mar = c(5.1,4.1,4.1,2.1)) 	
 	}
 	
