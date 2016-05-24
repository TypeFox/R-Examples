



##### Functions to identify sensible numbers of bins - help file on desktop  with data-frame setup
#####

# =============================================================================
# =============================================================================
convergeIPM<-function(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
		discreteTrans = 1, integrateType = "midpoint", correction = "none", 
		preCensus = TRUE, tol=1e-4,
		binIncrease=5, chosenBin=1, response="lambda"){
	
	if (response=="lambda") 
		rc <- .convergeLambda(growObj=growObj, survObj=survObj, fecObj=fecObj, 
				nBigMatrix=nBigMatrix, minSize=minSize, maxSize=maxSize, 
				discreteTrans =discreteTrans, integrateType = integrateType, 
				correction = correction, preCensus = preCensus,
				tol=tol,binIncrease=binIncrease)
			
	if (response=="R0")	
		rc <- .convergeR0(growObj=growObj, survObj=survObj, fecObj=fecObj, 
				nBigMatrix=nBigMatrix, minSize=minSize, maxSize=maxSize, 
				discreteTrans =discreteTrans, integrateType = integrateType, 
				correction = correction, preCensus = preCensus,
				tol=tol,binIncrease=binIncrease)
			
	if (response=="lifeExpect")	
		rc <- .convergeLifeExpectancy(growObj=growObj, survObj=survObj, 
			nBigMatrix=nBigMatrix, minSize=minSize, maxSize=maxSize, 
			discreteTrans =discreteTrans, integrateType = integrateType, 
			correction = correction, 
			tol=tol,binIncrease=binIncrease, chosenBin=chosenBin)


	return(rc)

}
	



### FUNCTIONS FOR EXTRACTING INTEGRATED DEMOGRAPHIC MEASURES ########################
### i.e. life-expectancy and passage time

# =============================================================================
# =============================================================================
#Generic for mean life expectancy
#parameters - an IPM
# returns - the life expectancy for every starting size. 
meanLifeExpect <- function(IPMmatrix){
	nBigMatrix <- length(IPMmatrix@.Data[1,]) #this nBigMatrix actually contains discrete, env, etc
	#tmp <-  ginv(diag(IPMmatrix@nEnvClass*nBigMatrix)-IPMmatrix)
	tmp <-  ginv(diag(nBigMatrix)-IPMmatrix)
	lifeExpect <- colSums(tmp)
	return(lifeExpect)
}


# =============================================================================
# =============================================================================
#Generic for variance life expectancy (p119 Caswell)
#parameters - an IPM
# returns - the variance in life expectancy for every starting size. 

varLifeExpect <- function(IPMmatrix){
	nBigMatrix <- length(IPMmatrix@.Data[1,])
	#tmp <-  ginv(diag(IPMmatrix@nEnvClass*nBigMatrix)-IPMmatrix)
	tmp <-  ginv(diag(nBigMatrix)-IPMmatrix)
	#varLifeExpect <- (2*diag(tmp)-diag(length(IPMmatrix[,1])))%*%tmp-(tmp*tmp)
	#varLifeExpect <- colSums(varLifeExpect)
	varLifeExpect <- colSums(2*(tmp%*%tmp)-tmp)-colSums(tmp)*colSums(tmp)                  
	return(varLifeExpect)
}




# =============================================================================
# =============================================================================
#Generic for survivorship
#parameters - IPMmatrix - an IPM
#           - size1 - a size at age 1
#           - maxAge - a maxAge
# returns - a list including the survivorship up to the max age,
#                      this broken down by stage,
#                       and mortality over age 


survivorship <- function(IPMmatrix, loc, maxAge=300){
	nBigMatrix <- length(IPMmatrix@.Data[1,])
	#n <- IPMmatrix@nEnvClass*nBigMatrix
	n <- nBigMatrix
	A1 <- tmp <-  IPMmatrix
	stage.agesurv <- matrix(NA,n,maxAge)
	surv.curv <- rep (NA,maxAge)
	
	#identify the starting size you want to track - removed - specify bin directly
	#loc <- which(abs(size1-IPMmatrix@meshpoints)==min(abs(size1-IPMmatrix@meshpoints)),arr.ind=T)[1]
	popvec <- matrix(0,n,1)
	popvec[floor(loc),1] <- 1
	
	for (a in 1:maxAge) {
		surv.curv[a]<-sum(A1[,loc]); 
		stage.agesurv[c(1:n),a]<-A1[,]%*%popvec
		A1<-A1%*%tmp
	}
	
	mortality <- -log(surv.curv[2:length(surv.curv)]/surv.curv[1:(length(surv.curv)-1)])
	
	return(list(surv.curv=surv.curv,stage.agesurv=stage.agesurv, mortality = mortality))
}

# =============================================================================
# =============================================================================
#Generic for first passage time 
#parameters - an IPM
#           - a size for which passage time is required            
# returns - the passage time to this size from each of the sizes in the IPM 

passageTime <- function(chosenSize,IPMmatrix){	
	loc <- which(abs(chosenSize-IPMmatrix@meshpoints) ==
					min(abs(chosenSize - IPMmatrix@meshpoints)),arr.ind=TRUE)[1]
	matrix.dim <- length(IPMmatrix[1,])
	
	Tprime <- IPMmatrix
	Tprime[,loc] <- 0
	
	Mprime <- 1-colSums(IPMmatrix)
	Mprime[loc]<-0
	Mprime <- rbind(Mprime,rep(0,matrix.dim))
	Mprime[2,loc] <- 1
	
	Bprime <- Mprime%*% ginv(diag(matrix.dim)-Tprime)
	#print(round(Bprime[2,],2))
	#print(sum(Bprime[2,]<1e-6))
	Bprime[2,][Bprime[2,]==0] <- 1
	
	diagBprime <- diag(Bprime[2,])
	#plot(IPMmatrix@meshpoints,diag(diagBprime),type="l",log="y")
	#abline(v=chosenSize)
	Tc <- diagBprime%*%Tprime%*%ginv(diagBprime)
	eta1 <- ginv(diag(matrix.dim)-Tc)             
	
	time.to.absorb <- colSums(eta1)
	time.to.absorb[loc:length(time.to.absorb)] <- 0
	return(time.to.absorb)
}

# =============================================================================
# =============================================================================
#Generic for first variance first passage time (not sure!!!)
#parameters - an IPM
#           - a size for which passage time is required            
# returns - the variance passage time to this size from each of the sizes in the IPM 
varPassageTime <- function(chosenSize,IPMmatrix){	
	loc <- which(abs(chosenSize-IPMmatrix@meshpoints)==min(abs(chosenSize-IPMmatrix@meshpoints)),arr.ind=TRUE)
	matrix.dim <- length(IPMmatrix[1,])
	
	Tprime <- IPMmatrix
	Tprime[,loc] <- 0
	
	Mprime <- 1-colSums(IPMmatrix)
	Mprime[loc]<-0
	Mprime <- rbind(Mprime,rep(0,matrix.dim))
	Mprime[2,loc] <- 1
	
	Bprime <- Mprime%*% solve(diag(matrix.dim)-Tprime)
	
	Tc <- diag(Bprime[2,])%*%Tprime%*%ginv(diag(Bprime[2,]))
	eta1 <- ginv(diag(matrix.dim)-Tc)             
	
	vartimeAbsorb <- colSums(2*(eta1%*%eta1)-eta1)-colSums(eta1)*colSums(eta1)                  
	
	return(vartimeAbsorb)
}

# =============================================================================
# =============================================================================
##Function to estimate Stochastic Passage Time
stochPassageTime <- function(chosenSize,IPMmatrix,envMatrix){
	#get the right index for the size you want
	loc <- which(abs(chosenSize-IPMmatrix@meshpoints)==min(abs(chosenSize-
									IPMmatrix@meshpoints)),arr.ind=TRUE)
	#expand out to find that size in every env
	#locs.all <- loc*c(1:IPMmatrix@nEnvClass)
	locs.all <- loc+((IPMmatrix@nBigMatrix)*(0:(envMatrix@nEnvClass-1)))
	matrix.dim <- length(IPMmatrix[1,])
	
	Tprime <- IPMmatrix
	Tprime[,locs.all] <- 0
	
	dhat <- 1-colSums(IPMmatrix)
	dhat[locs.all]<-0
	dhat <- rbind(dhat,rep(0,matrix.dim))
	dhat[2,locs.all] <- 1
	
	bhat <- dhat%*% solve(diag(matrix.dim)-Tprime)
	
	Mc <- diag(bhat[2,])%*%Tprime%*%solve(diag(bhat[2,]))
	eta1 <- solve(diag(matrix.dim)-Mc)             
	
	time.to.absorb <-colSums(eta1)
	
	return(time.to.absorb)
}



### Short term changes / matrix iteration functions ##############################################
# =============================================================================
# =============================================================================
## TO DO - ADJUST TO ALLOW DISCRETE CLASSES ##

## Function to predict distribution in x time-steps given starting
## distribution and IPM (with a single covariate)
#
#
#  Parameters - startingSizes - vector of starting sizes
#             - IPM the IPM (Pmatrix if only intrested in grow surv; Pmatrix + Fmatrix otherwise)
#             - n.time.steps - number of time steps
#             - startingEnv - vector of starting env, same length as startingSizes, or length=1
#
# Returns - a list including starting numbers in each IPM size class (n.new.dist0) and
#                            final numbers in each IPM size class (n.new.dist)
#
#
predictFutureDistribution <- function(startingSizes,IPM, n.time.steps, startingEnv=1) {
	
	# turn starting sizes into the resolution of the IPM bins
	breakpoints <- c(IPM@meshpoints-(IPM@meshpoints[2]-IPM@meshpoints[1]),
			IPM@meshpoints[length(IPM@meshpoints)]+(IPM@meshpoints[2]-IPM@meshpoints[1]))
	
	# setup slightly different for coompound or non compound dists
	if (IPM@nEnvClass>1) {
		if (length(startingEnv)==1) startingEnv <- rep(startingEnv, length(startingSizes))
		compound <- TRUE
		env.index <- IPM@env.index
		n.new.dist <- rep(0,length(IPM[1,]))
		for (ev in 1:IPM@nEnvClass) { 			
			index.new.dist <- findInterval(startingSizes[startingEnv==ev],breakpoints,all.inside=TRUE)
			loc.sizes <- table(index.new.dist); 
			n.new.dist[ev==IPM@env.index][as.numeric(names(loc.sizes))] <- loc.sizes
		}
		n.new.dist0 <- n.new.dist
	} else {
		compound <- FALSE
		index.new.dist <- findInterval(startingSizes,breakpoints,all.inside=TRUE)
		loc.sizes <- table(index.new.dist); 
		env.index <- rep(1,length(IPM@meshpoints))
		n.new.dist <- rep(0,length(IPM@meshpoints))
		n.new.dist[as.numeric(names(loc.sizes))] <- loc.sizes
		n.new.dist0 <- n.new.dist
	}
	
	for (t in 1:n.time.steps) n.new.dist <- IPM@.Data%*%n.new.dist
	
	plot(IPM@meshpoints,n.new.dist0[env.index==1],type="l",xlab="size",
			ylab="n in each size class", ylim=range(c(n.new.dist0,n.new.dist)))
	points(IPM@meshpoints,n.new.dist[env.index==1],type="l",col=2)
	if (compound) {
		for (j in 1:max(env.index)) {
			points(IPM@meshpoints,n.new.dist0[env.index==j],type="l",col=1,lty=j)
			points(IPM@meshpoints,n.new.dist[env.index==j],type="l",col=2,lty=j)
		}
		
	}
	legend("topright",legend=c("current","future"),col=1:2,lty=1,bty="n")
	
	
	return(list(n.new.dist0=n.new.dist0,n.new.dist=n.new.dist))     
}

# =============================================================================
# =============================================================================
## TO DO - ADJUST TO ALLOW DISCRETE CLASSES ##

## Function to see how long it takes to get from a starting distribution to an final size
##
#  Parameters - startingSizes - vector of starting sizes (in any order)
#             - IPM the IPM (just a Pmatrix)
#             - endSize - the end size
#             - startingEnv - vector of starting env, same length as startingSizes, or length=1
#             - maxT - the max number of time-steps tested
#             - propReach - the proportion of the starting pop that have to be > than the endSize for it to count
#                  (plots and returned values of survivorship from preliminary runs will give a notion of how low this has to be)
#
# Returns - a list containing: ts.dist - the time-series of size distribution
#                              time.reach - the time for n.reach to be at sizes > endSize
#                              survivorship - survivorship over the course of the time elapsed for that pop
#                              
#
# THE PROBLEM WITH THIS FUNCTION IS THAT EITHER 1) YOU MAKE IT IN ONE TIME-STEP; OR 2) EVERYONE IS DEAD SO TUNING
# propReach BECOMES THE KEY - AND THE EXACT VALUE TO PROVIDE VALUES 1 < x < maxT CAN BE LUDICRIOUSLY SENSITIVE
#
timeToSize <- function(startingSizes,IPM,endSize, startingEnv=1, maxT=100, propReach=0.01) {
	
	cutoff <- which(IPM@meshpoints>endSize,arr.ind=TRUE)[1]
	n.reach <- propReach*length(startingSizes)
	
	# setup slightly different for coompound or non compound dists
	if (IPM@nEnvClass>1) {
		#if startingEnv is not a vector, assume all start in startingEnv
		if (length(startingEnv)==1) startingEnv <- rep(startingEnv, length(startingSizes))
		compound <- TRUE
		env.index <- IPM@env.index
		n.new.dist <- rep(0,length(IPM[1,]))
		for (ev in 1:IPM@nEnvClass) { 
			index.new.dist <- findInterval(startingSizes[startingEnv==ev],IPM@meshpoints)+1
			index.new.dist[index.new.dist>length(IPM@meshpoints)] <- length(IPM@meshpoints)
			loc.sizes <- table(index.new.dist); 
			n.new.dist[ev==IPM@env.index][as.numeric(names(loc.sizes))] <- loc.sizes
		}
		n.new.dist0 <- n.new.dist
	} else {
		compound <- FALSE
		index.new.dist <- findInterval(startingSizes,IPM@meshpoints)+1
		index.new.dist[index.new.dist>length(IPM@meshpoints)] <- length(IPM@meshpoints)
		loc.sizes <- table(index.new.dist); 
		env.index <- rep(1,length(IPM@meshpoints))
		n.new.dist <- rep(0,length(IPM@meshpoints))
		n.new.dist[as.numeric(names(loc.sizes))] <- loc.sizes
		n.new.dist0 <- n.new.dist
	}
	
	ts.dist <- matrix(NA,length(n.new.dist),maxT)
	
	survivorship <- rep(NA,maxT)
	for (t in 1:maxT) { 
		#print(t)
		n.new.dist <- IPM@.Data%*%n.new.dist
		#plot(n.new.dist)
		ts.dist[,t] <- n.new.dist
		
		if (!compound) {
			tot <- sum(n.new.dist[cutoff:length(IPM@meshpoints)])
			survivorship[t] <- sum(n.new.dist)/length(startingSizes)
		} else {
			tot <-sumN <- 0
			for (ev in 1:IPM@nEnvClass) {
				tot <- tot+sum(n.new.dist[env.index==ev][cutoff:length(IPM@meshpoints)])
				sumN <- sumN + sum(n.new.dist[env.index==ev])
			}
			survivorship[t] <- sumN/length(startingSizes)
		}
		
		
		if (tot>n.reach){
			time.reach <- t
			break()
		}
		
	}
	
	if (t==maxT) time.reach <- maxT
	
	par(mfrow=c(1,3),bty="l")
	plot(IPM@meshpoints,n.new.dist0[env.index==1],type="l",xlab="size", ylab="n", ylim=range(c(n.new.dist0+10,n.new.dist)))
	points(IPM@meshpoints,n.new.dist[env.index==1],type="l",col=2)
	legend("topleft",legend=c("starting distribution", "final distribution"),col=c(1,2),lty=1,bty="n")
	abline(v=IPM@meshpoints[cutoff],lty=3)
	
	plot(survivorship[1:t], xlab="Time", ylab="survivorship", type="l")
	
	if (time.reach>5) { 
		image(as.numeric(IPM@meshpoints),1:time.reach,log(ts.dist),ylab="Time steps", xlab="Size classes", main="numbers in size classes over time")
		contour(as.numeric(IPM@meshpoints),1:time.reach,log(ts.dist),add=TRUE,levels=exp(seq(0,max(log(ts.dist)),length=10)))
	}
	print(paste("Time to reach:",time.reach))
	
	return(list(ts.dist=ts.dist, time.reach=time.reach, survivorship=survivorship))     
}







# =============================================================================
# =============================================================================
# Calculate R0 
#
# Parameters - Fmatrix, Pmatrix
#
# Returns R0
R0Calc<-function(Pmatrix, Fmatrix){
	Imatrix <- matrix(0, length(Pmatrix[1,]), length(Pmatrix[1,])); 
	diag(Imatrix) <- 1
	Nmatrix <- ginv(Imatrix - Pmatrix);
	Rmatrix <- Fmatrix %*% Nmatrix
	ave.R0 <- Re(eigen(Rmatrix)$values[1])
	return(ave.R0)
}

# =============================================================================
# =============================================================================
# Calculate lambda and stable stage dist
# for constant env for really huge matrices
# using Matrix package for numerical efficiency
#
# Parameters - Amat - an IPM object
#            - tol - tolerance (i.e. precision required)
#
# Returns list containing lambda and stableStage
#
#ROB WILL MODIFY THIS CODE TO INCLUDE REPRODUCTIVE VALUES
largeMatrixCalc <- function(Pmatrix, Fmatrix, tol = 1.e-8){
	A2 <- Matrix(Pmatrix + Fmatrix);
	nt <- Matrix(1,length(Pmatrix[1,]), 1);
	nt1 <- nt; 
	
	h1 <- diff(Pmatrix@meshpoints)[1]
	
	qmax <- 1000;
	lam <- 1; 
	while(qmax > tol) {
		nt1 <- A2 %*% nt;
		qmax <- sum(abs((nt1 - lam * nt)@x));  
		lam <- sum(nt1@x); 
		nt@x <- (nt1@x) / lam; #slight cheat  
		#cat(lam,qmax,"\n");
	} 
	nt <- matrix(nt@x, length(Pmatrix[1,]), 1); 
	stableDist <- nt / (h1 * sum(nt)); #normalize so that integral=1
	lamStable <- lam; 
	
	# Check works   
	qmax <- sum(abs(lam * nt - (Pmatrix + Fmatrix) %*% nt))
	cat("Convergence: ", qmax, " should be less than ", tol, "\n")
	
	
	return(list(lam = lam, stableDist = stableDist, h1 = h1)) 
	
}

# =============================================================================
# =============================================================================
## Sensitivity of parameters - works for an IPM built out of
## growth, survival, discreteTrans, fecundity and clonality objects.
##
##
sensParams <- function (growObj, survObj, fecObj=NULL, clonalObj=NULL,
		nBigMatrix, minSize, maxSize,
		chosenCov = data.frame(covariate = 1), discreteTrans = 1,
		integrateType = "midpoint", correction = "none", 
		preCensusFec = TRUE, postCensusSurvObjFec = NULL, postCensusGrowObjFec = NULL,  
		preCensusClonal = TRUE, postCensusSurvObjClonal = NULL, postCensusGrowObjClonal = NULL,  
		delta = 1e-04, response="lambda", chosenBin=1) {
	
	if (response!="lambda" & response!="R0" & response !="lifeExpect")
		stop("response must be one of lambda or R0 or lifeExpect")
	
	nmes <- elam <- slam <- c()
	
	# get the base
	Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
			chosenCov = chosenCov, maxSize = maxSize, growObj = growObj,
			survObj = survObj, discreteTrans = discreteTrans, integrateType = integrateType,
			correction = correction)
	if (!is.null(fecObj)) {
		Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				chosenCov = chosenCov, maxSize = maxSize, fecObj = fecObj,
				integrateType = integrateType, correction = correction,
				preCensus = preCensusFec, survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
	} else {Fmatrix <- Pmatrix*0 }
	if (!is.null(clonalObj)) {
		Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
				integrateType = integrateType, correction = correction,
				preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
				growObj = postCensusGrowObjClonal)
	} else {Cmatrix <- Pmatrix*0 }
	
	IPM <- Pmatrix + Fmatrix + Cmatrix
	
	if (response=="lambda") rc1 <- Re(eigen(IPM)$value[1])
	if (response=="R0") rc1 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
	if (response=="lifeExpect") rc1 <- meanLifeExpect(Pmatrix)[chosenBin]
	
	# 1. survival
	for (j in 1:length(survObj@fit$coeff)) {
		survObj@fit$coefficients[j] <- survObj@fit$coefficients[j] * (1 + delta)
		Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
				minSize = minSize, maxSize = maxSize, growObj = growObj,
				survObj = survObj, discreteTrans = discreteTrans,
				chosenCov = chosenCov, integrateType = integrateType,
				correction = correction)
		if (!is.null(fecObj)) {
			Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
					minSize = minSize, maxSize = maxSize, fecObj = fecObj,
					integrateType = integrateType, correction = correction,
					chosenCov = chosenCov, preCensus = preCensusFec, 
					survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
		} 
		if (!is.null(clonalObj)) {
			Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
					integrateType = integrateType, correction = correction,
					preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
					growObj = postCensusGrowObjClonal)
		} 
		
		IPM <- Pmatrix + Fmatrix + Cmatrix
		
		if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
		if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
		if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
		
		survObj@fit$coefficients[j] <- survObj@fit$coefficients[j]/(1 + delta)
		
		slam <- c(slam, (rc2 - rc1)/((as.numeric(survObj@fit$coefficients[j]))* delta))
		elam <- c(elam, (rc2 - rc1)/(rc1 *delta))
		nmes <- c(nmes, as.character(paste("survival:",names(survObj@fit$coeff)[j])))
	}
	
	# 2 growth
	for (j in 1:length(growObj@fit$coeff)) {
		growObj@fit$coefficients[j] <- growObj@fit$coefficients[j] * (1 + delta)
		Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
				minSize = minSize, maxSize = maxSize, growObj = growObj,
				chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
				integrateType = integrateType, correction = correction)
		if (!is.null(fecObj)) {
			Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
					minSize = minSize, maxSize = maxSize, fecObj = fecObj,
					chosenCov = chosenCov, integrateType = integrateType,
					correction = correction, preCensus = preCensusFec, 
					survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
		}
		if (!is.null(clonalObj)) {
			Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
					integrateType = integrateType, correction = correction,
					preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
					growObj = postCensusGrowObjClonal)
		}
		
		IPM <- Pmatrix + Fmatrix + Cmatrix
		
		if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
		if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
		if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
		
		growObj@fit$coefficients[j] <- growObj@fit$coefficients[j]/(1 + delta)
		
		slam <- c(slam, (rc2 - rc1)/(as.numeric(growObj@fit$coefficients[j]) * delta))
		elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
		nmes <- c(nmes, as.character(paste("growth:",names(growObj@fit$coeff)[j])))
	}
	
	growObj@sd <- growObj@sd * (1 + delta)
	Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
			maxSize = maxSize, growObj = growObj, survObj = survObj,
			chosenCov = chosenCov, discreteTrans = discreteTrans,
			integrateType = integrateType, correction = correction)
	if (!is.null(fecObj)) {
		Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
				chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
				survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
	}
	if (!is.null(clonalObj)) {
		Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
				integrateType = integrateType, correction = correction,
				preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
				growObj = postCensusGrowObjClonal)
	}
	IPM <- Pmatrix + Fmatrix + Cmatrix
	
	if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
	if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
	if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
	
	growObj@sd <- growObj@sd / (1 + delta)
	
	slam <- c(slam,(rc2 - rc1)/(growObj@sd * delta))
	elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
	nmes <- c(nmes, "growth: sd")
	
	# 3. DiscreteTrans
	if (class(discreteTrans)=="discreteTrans") {
		for (j in 1:(ncol(discreteTrans@discreteTrans)-1)) {
			for (i in 1:(nrow(discreteTrans@discreteTrans)-1)) {
				discreteTrans@discreteTrans[i,j]<-discreteTrans@discreteTrans[i,j] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						maxSize = maxSize, growObj = growObj, survObj = survObj,
						chosenCov = chosenCov, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
							chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				discreteTrans@discreteTrans[i,j]<-discreteTrans@discreteTrans[i,j] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(discreteTrans@discreteTrans[i,j] * delta))
				elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes, as.character(paste("discrete:",dimnames(discreteTrans@discreteTrans)[[2]][j],"to",dimnames(discreteTrans@discreteTrans)[[1]][i])))
			}
		}
		#if there is more than 2 discrete stages (beyond "continuous" "dead" and one discrete stage)
        #then moveToDiscrete tells you how many of surviving continuous individuals are going into 
		#discrete classes, but how they distributed also; which is the last column in discreteTrans 
		if (nrow(discreteTrans@discreteTrans)>3) {
			for (i in 1:(nrow(discreteTrans@discreteTrans)-2)) {
				discreteTrans@discreteTrans[i,"continuous"]<-discreteTrans@discreteTrans[i,"continuous"] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						maxSize = maxSize, growObj = growObj, survObj = survObj,
						chosenCov = chosenCov, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
							chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				discreteTrans@discreteTrans[i,"continuous"]<-discreteTrans@discreteTrans[i,"continuous"] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(discreteTrans@discreteTrans[i,"continuous"] * delta))
				elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes, as.character(paste("discrete: Continuous to",dimnames(discreteTrans@discreteTrans)[[1]][i])))
			}
		}
		
		for (j in 1:length(discreteTrans@meanToCont)) {
			discreteTrans@meanToCont[1,j]<-discreteTrans@meanToCont[1,j] * (1 + delta)
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, growObj = growObj, survObj = survObj,
					chosenCov = chosenCov, discreteTrans = discreteTrans,
					integrateType = integrateType, correction = correction)
			if (!is.null(fecObj)) {
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
						chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
			}
			if (!is.null(clonalObj)) {
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
			}
			IPM <- Pmatrix + Fmatrix + Cmatrix
			
			if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
			if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
			if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
			
			discreteTrans@meanToCont[1,j]<-discreteTrans@meanToCont[1,j] / (1 + delta)
			
			slam <- c(slam,(rc2 - rc1)/(discreteTrans@meanToCont[1,j] * delta))
			elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
			nmes <- c(nmes, as.character(paste("discrete: meanToCont",dimnames(discreteTrans@meanToCont)[[2]][j])))
		}
		
		for (j in 1:length(discreteTrans@sdToCont)) {
			discreteTrans@sdToCont[1,j]<-discreteTrans@sdToCont[1,j] * (1 + delta)
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, growObj = growObj, survObj = survObj,
					chosenCov = chosenCov, discreteTrans = discreteTrans,
					integrateType = integrateType, correction = correction)
			if (!is.null(fecObj)) {
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
						chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
			}
			if (!is.null(clonalObj)) {
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
			}
			IPM <- Pmatrix + Fmatrix + Cmatrix
			
			if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
			if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
			if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
			
			discreteTrans@sdToCont[1,j]<-discreteTrans@sdToCont[1,j] / (1 + delta)
			
			slam <- c(slam,(rc2 - rc1)/(discreteTrans@sdToCont[1,j] * delta))
			elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
			nmes <- c(nmes, as.character(paste("discrete: sdToCont",dimnames(discreteTrans@sdToCont)[[2]][j])))
		}
		
		for (j in 1:length(discreteTrans@moveToDiscrete$coef)) {
			discreteTrans@moveToDiscrete$coefficients[j]<-discreteTrans@moveToDiscrete$coefficients[j] * (1 + delta)
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, growObj = growObj, survObj = survObj,
					chosenCov = chosenCov, discreteTrans = discreteTrans,
					integrateType = integrateType, correction = correction)
			if (!is.null(fecObj)) {
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						maxSize = maxSize, fecObj = fecObj, integrateType = integrateType,
						chosenCov = chosenCov, correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
			}
			if (!is.null(clonalObj)) {
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
			}
			IPM <- Pmatrix + Fmatrix + Cmatrix
			
			if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
			if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
			if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
			
			discreteTrans@moveToDiscrete$coefficients[j]<-discreteTrans@moveToDiscrete$coefficients[j] / (1 + delta)
			
			slam <- c(slam,(rc2 - rc1)/(as.numeric(discreteTrans@moveToDiscrete$coefficients[j]) * delta))
			elam <- c(elam, (rc2 - rc1)/(rc1 * delta))
			nmes <- c(nmes, as.character(paste("discrete: moveToDiscrete",names(discreteTrans@moveToDiscrete$coefficients)[j])))
		}
	}
	
	# 4. Fecundity
	if (!is.null(fecObj)) {
		for (i in 1:length(fecObj@fitFec)) {
			for (j in 1:length(fecObj@fitFec[[i]]$coefficients)) {
				fecObj@fitFec[[i]]$coefficients[j] <- fecObj@fitFec[[i]]$coefficients[j] *
						(1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						chosenCov = chosenCov, integrateType = integrateType,
						correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				fecObj@fitFec[[i]]$coefficients[j] <- fecObj@fitFec[[i]]$coefficients[j]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/((as.numeric(fecObj@fitFec[[i]]$coefficients[j]) * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("fecundity: func", i, names(fecObj@fitFec[[i]]$coefficients)[j]))
			}
		}
		
		chs <- which(!is.na(as.numeric(fecObj@fecConstants)), arr.ind = TRUE)
		if (length(chs) > 0) {
			for (j in 1:length(chs)) {
				fecObj@fecConstants[1,chs[j]] <- fecObj@fecConstants[1,chs[j]] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						chosenCov = chosenCov, integrateType = integrateType,
						correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				fecObj@fecConstants[1, chs[j]] <- fecObj@fecConstants[1,chs[j]]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(fecObj@fecConstants[1,chs[j]] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("fecundity: constant",names(fecObj@fecConstants)[chs[j]]))
			}
		}
		
		if (max(fecObj@offspringSplitter)<1) {
			for (j in which(fecObj@offspringSplitter>0)) {
				fecObj@offspringSplitter[j] <- fecObj@offspringSplitter[j] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						chosenCov = chosenCov, integrateType = integrateType,
						correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				fecObj@offspringSplitter[j] <- fecObj@offspringSplitter[j] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(fecObj@offspringSplitter[j] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("fecundity: offspringSplitter",names(fecObj@offspringSplitter[j])))
			}
		}
		
		chs <- which(!is.na(as.numeric(fecObj@fecByDiscrete)), arr.ind = TRUE)
		if (length(chs) > 0) {
			for (j in 1:length(chs)) {
				fecObj@fecByDiscrete[1,chs[j]] <- fecObj@fecByDiscrete[1,chs[j]] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						chosenCov = chosenCov, integrateType = integrateType,
						correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				fecObj@fecByDiscrete[1,chs[j]] <- fecObj@fecByDiscrete[1,chs[j]] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(fecObj@fecByDiscrete[1,chs[j]] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("fecundity: fecByDiscrete",names(fecObj@fecByDiscrete)[chs[j]]))
			}
		}
		
		if (class(fecObj@offspringRel)=="lm") {
			for (j in 1:length(fecObj@offspringRel$coeff)) {
				fecObj@offspringRel$coefficients[j] <- fecObj@offspringRel$coefficients[j] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						integrateType = integrateType, correction = correction,
						chosenCov = chosenCov, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				if (!is.null(clonalObj)) {
					Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
							chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
							integrateType = integrateType, correction = correction,
							preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
							growObj = postCensusGrowObjClonal)
				}
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				fecObj@offspringRel$coefficients[j] <- fecObj@offspringRel$coefficients[j]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(fecObj@offspringRel$coefficients[j]) *delta))
				elam <- c(elam,(rc2 - rc1)/(rc1 *delta))
				nmes <- c(nmes, as.character(paste("fecundity: offspring rel ",names(fecObj@offspringRel$coeff)[j])))
			}
			
			fecObj@sdOffspringSize <- fecObj@sdOffspringSize * (1 + delta)
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, growObj = growObj, chosenCov = chosenCov,
					survObj = survObj, discreteTrans = discreteTrans, integrateType = integrateType,
					correction = correction)
			Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, fecObj = fecObj, chosenCov = chosenCov,
					integrateType = integrateType, correction = correction,
					preCensus = preCensusFec, 
					survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
			if (!is.null(clonalObj)) {
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
			}
			
			IPM <- Pmatrix + Fmatrix + Cmatrix
			
			if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
			if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
			if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
			
			fecObj@sdOffspringSize <- fecObj@sdOffspringSize/(1 + delta)
			
			slam <- c(slam,(rc2 - rc1)/(fecObj@sdOffspringSize * delta))
			elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
			nmes <- c(nmes, "fecundity: sd offspring size")
		}
	}
	
# 5. Clonality
	if (!is.null(clonalObj)) {
		for (i in 1:length(clonalObj@fitFec)) {
			for (j in 1:length(clonalObj@fitFec[[i]]$coefficients)) {
				clonalObj@fitFec[[i]]$coefficients[j] <- clonalObj@fitFec[[i]]$coefficients[j] *
						(1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
							minSize = minSize, maxSize = maxSize, fecObj = fecObj,
							chosenCov = chosenCov, integrateType = integrateType,
							correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				clonalObj@fitFec[[i]]$coefficients[j] <- clonalObj@fitFec[[i]]$coefficients[j]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/((as.numeric(clonalObj@fitFec[[i]]$coefficients[j]) * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("clonality: func", i, names(clonalObj@fitFec[[i]]$coefficients)[j]))
			}
		}
		
		chs <- which(!is.na(as.numeric(clonalObj@fecConstants)), arr.ind = TRUE)
		if (length(chs) > 0) {
			for (j in 1:length(chs)) {
				clonalObj@fecConstants[1,chs[j]] <- clonalObj@fecConstants[1,chs[j]] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
							minSize = minSize, maxSize = maxSize, fecObj = fecObj,
							chosenCov = chosenCov, integrateType = integrateType,
							correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				clonalObj@fecConstants[1, chs[j]] <- clonalObj@fecConstants[1,chs[j]]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(clonalObj@fecConstants[1,chs[j]] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("clonality: constant",names(clonalObj@fecConstants)[chs[j]]))
			}
		}
		
		if (max(clonalObj@offspringSplitter)<1) {
			for (j in which(clonalObj@offspringSplitter>0)) {
				clonalObj@offspringSplitter[j] <- clonalObj@offspringSplitter[j] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
							minSize = minSize, maxSize = maxSize, fecObj = fecObj,
							chosenCov = chosenCov, integrateType = integrateType,
							correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				clonalObj@offspringSplitter[j] <- clonalObj@offspringSplitter[j] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(clonalObj@offspringSplitter[j] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("clonality: offspringSplitter",names(clonalObj@offspringSplitter[j])))
			}
		}
		
		chs <- which(!is.na(as.numeric(clonalObj@fecByDiscrete)), arr.ind = TRUE)
		if (length(chs) > 0) {
			for (j in 1:length(chs)) {
				clonalObj@fecByDiscrete[1,chs[j]] <- clonalObj@fecByDiscrete[1,chs[j]] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
							minSize = minSize, maxSize = maxSize, fecObj = fecObj,
							chosenCov = chosenCov, integrateType = integrateType,
							correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				clonalObj@fecByDiscrete[1,chs[j]] <- clonalObj@fecByDiscrete[1,chs[j]] / (1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(clonalObj@fecByDiscrete[1,chs[j]] * delta)))
				elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
				nmes <- c(nmes,paste("clonality: fecByDiscrete",names(clonalObj@fecByDiscrete)[chs[j]]))
			}
		}
		
		if (class(clonalObj@offspringRel)=="lm") {
			for (j in 1:length(clonalObj@offspringRel$coeff)) {
				clonalObj@offspringRel$coefficients[j] <- clonalObj@offspringRel$coefficients[j] * (1 + delta)
				Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, growObj = growObj,
						chosenCov = chosenCov, survObj = survObj, discreteTrans = discreteTrans,
						integrateType = integrateType, correction = correction)
				if (!is.null(fecObj)) {
					Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
							minSize = minSize, maxSize = maxSize, fecObj = fecObj,
							chosenCov = chosenCov, integrateType = integrateType,
							correction = correction, preCensus = preCensusFec, 
							survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
				}
				Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
						chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
						integrateType = integrateType, correction = correction,
						preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
						growObj = postCensusGrowObjClonal)
				
				IPM <- Pmatrix + Fmatrix + Cmatrix
				
				if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
				if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
				if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
				
				clonalObj@offspringRel$coefficients[j] <- clonalObj@offspringRel$coefficients[j]/(1 + delta)
				
				slam <- c(slam,(rc2 - rc1)/(as.numeric(clonalObj@offspringRel$coefficients[j]) *delta))
				elam <- c(elam,(rc2 - rc1)/(rc1 *delta))
				nmes <- c(nmes, as.character(paste("clonality: offspring rel ",names(clonalObj@offspringRel$coeff)[j])))
			}
			
			clonalObj@sdOffspringSize <- clonalObj@sdOffspringSize * (1 + delta)
			Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					maxSize = maxSize, growObj = growObj, chosenCov = chosenCov,
					survObj = survObj, discreteTrans = discreteTrans, integrateType = integrateType,
					correction = correction)
			if (!is.null(fecObj)) {
				Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix,
						minSize = minSize, maxSize = maxSize, fecObj = fecObj,
						chosenCov = chosenCov, integrateType = integrateType,
						correction = correction, preCensus = preCensusFec, 
						survObj = postCensusSurvObjFec, growObj = postCensusGrowObjFec)
			}
			Cmatrix <- makeIPMCmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
					chosenCov = chosenCov, maxSize = maxSize, clonalObj = clonalObj,
					integrateType = integrateType, correction = correction,
					preCensus = preCensusClonal, survObj = postCensusSurvObjClonal, 
					growObj = postCensusGrowObjClonal)
			
			IPM <- Pmatrix + Fmatrix + Cmatrix
			
			if (response=="lambda") rc2 <- Re(eigen(IPM)$value[1])
			if (response=="R0") rc2 <- R0Calc(Pmatrix, Fmatrix+Cmatrix)
			if (response=="lifeExpect") rc2 <- meanLifeExpect(Pmatrix)[chosenBin]
			
			clonalObj@sdOffspringSize <- clonalObj@sdOffspringSize/(1 + delta)
			
			slam <- c(slam,(rc2 - rc1)/(clonalObj@sdOffspringSize * delta))
			elam <- c(elam,(rc2 - rc1)/(rc1 * delta))
			nmes <- c(nmes, "clonality: sd offspring size")
		}
	}
	
	names(slam) <- nmes
	names(elam) <- nmes
	
	return(list(sens = slam, elas = elam))
}


### FUNCTIONS FOR EXTRACTING STOCH GROWTH RATE ########################
# =============================================================================
# =============================================================================
# Generic approach to get stoch rate
# by sampling list IPM
#
# Parameters - listIPMmatrix - list of IPMs corresponding to different year types
#            - nRunIn - the burnin before establishing lambda_s
#            - tMax - the total time-horizon for getting lambda_s
#
# Returns lambda_s (no density dependence)

stochGrowthRateSampleList <- function(nRunIn,tMax,listIPMmatrix=NULL,
					listPmatrix=NULL, listFmatrix=NULL, seedList = NULL,
					densDep=FALSE){
	
	if (densDep & (is.null(listPmatrix) | is.null(listFmatrix))){
		stop("Require listPmatrix & listFmatrix for densDep=TRUE")
	}
		
	if (!densDep & is.null(listIPMmatrix)) {
		stop("Require listIPMmatrix for densDep=TRUE")		
	}		
	
	if (densDep) {
		nmatrices <- length(listPmatrix)
		nBigMatrix <- length(listPmatrix[[1]][,1]) 
	} else { 
		nmatrices <- length(listIPMmatrix)
		nBigMatrix <- length(listIPMmatrix[[1]][,1])
	}
	
	nt<-rep(1,nBigMatrix)
	Rt<-rep(NA,tMax)
	
	pEst <- 1
	
	for (t in 1:tMax) {
		year.type <- sample(1:nmatrices,size=1,replace=FALSE)
		if (densDep) { 
			nseeds <- sum(listFmatrix[[year.type]]%*%nt)
			pEst <- min(seedList[min(year.type+1,nmatrices)]/nseeds,1)
			nt1 <- (listPmatrix[[year.type]]+pEst*listFmatrix[[year.type]])%*% nt
			sum.nt1 <- sum(nt1)
			Rt[t] <- log(sum.nt1)-log(sum(nt))
			nt <- nt1
		} else {
			nt1<-listIPMmatrix[[year.type]] %*% nt	
			sum.nt1<-sum(nt1)
			Rt[t]<-log(sum.nt1)
			nt<-nt1/sum.nt1
		}		
	}
		res <- mean(Rt[nRunIn:tMax],na.rm=TRUE)
	return(res)
}

# =============================================================================
# =============================================================================
# Approach to get stoch rate
# with time-varying covariates
#
# Parameters - covariate - the covariates (temperature, etc)
#            - nRunIn - the burnin before establishing lambda_s
#            - tMax - the total time-horizon for getting lambda_s
#
# Returns lambda_s 


stochGrowthRateManyCov <- function(covariate,nRunIn,tMax,
		growthObj,survObj,fecObj,
		nBigMatrix,minSize,maxSize, nMicrosites,
		integrateType="midpoint",correction="none", 
		trackStruct=FALSE, plot=FALSE,...){	
	
	nt<-rep(1,nBigMatrix)
	Rt<-rep(NA,tMax)
	fecObj@fecConstants[is.na(fecObj@fecConstants)] <- 1 
	tmp.fecObj <- fecObj
	
	#print("fec.const start")
	#print(fecObj@fecConstants)
	
	#density dep in seedling establishment 
	if (sum(nMicrosites)>0) { dd <- TRUE; seeds <- 10000 } else { dd <- FALSE; p.est <- 1}
	if (trackStruct) rc <- matrix(NA,nBigMatrix,tMax)
	
	for (t in 1:tMax) {
		if (dd) p.est <- min(nMicrosites[min(t,length(nMicrosites))]/seeds,1) 	
		
		#if you don't do this, rownames return errors...
		covariatet <- covariate[t,]
		row.names(covariatet) <- NULL
		
		#but if you have only one column, then it can forget its a data-frame
		if (ncol(covariate)==1) { 
			covariatet <- data.frame(covariatet)
			colnames(covariatet) <- colnames(covariate)	
		}
		
		tpS <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = covariatet,
				growObj = growthObj, survObj = survObj,
				integrateType=integrateType, correction=correction)
		
		tpF <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = covariatet,
				fecObj = tmp.fecObj,
				integrateType=integrateType, correction=correction)
		
		#total seeds for next year 
		if (dd) seeds <- sum(tpF%*%nt)
		
		IPM.here <- p.est*tpF@.Data+tpS@.Data
		nt1<-IPM.here %*% nt	
		sum.nt1<-sum(nt1)
		
		if (!trackStruct){
			if (!dd) { 
				Rt[t]<-log(sum.nt1)
				nt<-nt1/sum.nt1
			} else {
				Rt[t]<-log(sum.nt1)-log(sum(nt))
				nt<-nt1	
			}} else {
			nt<-nt1	
			rc[,t] <- nt1
		}
		
		
		
	}
	
	if (trackStruct & plot){
		.plotResultsStochStruct(tVals=1:tMax, meshpoints=tpS@meshpoints,
				st=rc,covTest=covariate, nRunIn = nRunIn,  ...) 	
	}
	
	if (!trackStruct) {res <- mean(Rt[nRunIn:tMax],na.rm=TRUE); return(list(Rt=res))}
	if (trackStruct) return(list(rc=rc,IPM.here=tpS))
	
}






## FUNCTIONS FOR SENSITIVITY AND ELASTICITY ###################################

# =============================================================================
# =============================================================================
#parameters - an IPM (with full survival and fecundity complement)
# returns - the sensitivity of every transition for pop growth
sens<-function(A) {
	w<-Re(eigen(A)$vectors[,1]); 
	v<-Re(eigen(t(A))$vectors[,1]);
	vw<-sum(v*w);
	s<-outer(v,w)
	return(s/vw); 
}   

# =============================================================================
# =============================================================================
#parameters - an IPM (with full survival and fecundity complement)
# returns - the elasticity of every transition for pop growth
elas<-function(A) {
	s<-sens(A)
	lam<-Re(eigen(A)$values[1]);
	return((s*A)/lam);
}


# =============================================================================
# =============================================================================
## Function to get passage time FROM a particular size TO a range of sizes
## (i.e. size to age) when provided with a Pmatrix, a starting size, and a list
## of target sizes
#
# Parameters - Pmatrix
#            - startingSize
#            - targetSizes
#
# Returns - list containing vector of targets, vector of corresponding times, and the startingSize
#
sizeToAge <- function(Pmatrix,startingSize,targetSize) {
  
  #locate where the first size is in the meshpoints of Pmatrix
  diffv <- abs(startingSize-Pmatrix@meshpoints)
  start.index <- median(which(diffv==min(diffv),arr.ind=TRUE))
  timeInYears <- rep(NA,length(targetSize))
  
  #loop over to see where its going
  for (k in 1:length(targetSize)) {
    pTime <- passageTime(targetSize[k],Pmatrix)
    timeInYears[k] <- pTime[start.index]
  }
  
  return(list(timeInYears=timeInYears,targetSize=targetSize,startingSize=startingSize))
  
}

