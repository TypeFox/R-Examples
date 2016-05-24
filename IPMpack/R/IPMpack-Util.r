


# =============================================================================
# =============================================================================
## FUNCTION FOR MAKING PICTURE OF SURVIVAL DATA AND FITTED SURVIVAL OBJECT ####
## and PICTURE GROWTH DATA AND FITTED GROWTH OBJECT  
# ** note that growth may need redefining if new growth  objects are defined


# parameters - dataf - the dataframe
#            - survObj - a survival object
#            - ncuts - the number of cuts used for binning survival
# returns - 

picSurv <- function(dataf, survObj, ncuts = 20, makeTitle = "Survival", ...) { 
	
	#organize data and plot mean of ncut successive sizes, so trend more obvious
	os<-order(dataf$size); os.surv<-(dataf$surv)[os]; os.size<-(dataf$size)[os]; 
	psz<-tapply(os.size,as.numeric(cut(os.size,ncuts)),mean,na.rm=TRUE); #print(psz)
	ps<-tapply(os.surv,as.numeric(cut(os.size,ncuts)), mean, na.rm = TRUE);#print(ps)
	
	if (length(grep("covariate",names(survObj@fit$model)))==0) {  
		#plot data
		plot(as.numeric(psz),as.numeric(ps),pch=19,
				xlab="Size at t", ylab = "Survival to t+1", main = makeTitle, ...)
		#Plot fitted models
		points(dataf$size[order(dataf$size)],surv(dataf$size[order(dataf$size)],data.frame(covariate=1),survObj),type="l",col=2)
	} else {
		plot(as.numeric(psz),as.numeric(ps),
				type="n",pch=19,xlab="Size at t", ylab="Survival to t+1",main="Survival",...)
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		os.cov<-(dataf$covariate)[os]
		sizes <- dataf$size[!is.na(dataf$size)]; sizes <- sizes[order(sizes)]
		ud <- unique(dataf$covariate); ud <- ud[!is.na(ud)]
		for (k in 1:length(ud)) { 
			tp <- os.cov==ud[k]
			psz<-tapply(os.size[tp], as.numeric(cut(os.size[tp], ncuts)), mean, na.rm = TRUE); #print(psz)
			ps<-tapply(os.surv[tp],as.numeric(cut(os.size[tp],ncuts)),mean,na.rm=TRUE);#print(ps)
			points(as.numeric(psz),as.numeric(ps),pch=19,col=k)
			newd <- data.frame(size=sizes,size2=sizes^2,size3=sizes^3,
					covariate=rep(as.factor(ud[k]),length(sizes)))
			if(length(grep("expsize",survObj@fit$formula))==1)
				newd$expsize=exp(sizes)
			if(length(grep("logsize",survObj@fit$formula))==1)
				newd$logsize=log(sizes)
			if(length(grep("logsize2",survObj@fit$formula))==1)
				newd$logsize=(log(sizes))^2
			
			pred.surv <- predict(survObj@fit,newd,type="response")
			points(newd$size,pred.surv,type="l",col=k)
		}
	} 
	

}


# =============================================================================
# =============================================================================
# Function to make pic of growth object fit with data
#
# parameters - dataf - the dataframe
#            - growObj - a growth object
#
# returns - 
#

picGrow <- function(dataf, growObj, mainTitle = "Growth",...) {
	predVar <- attr(growObj@fit$terms,"predvars")[[2]]  #jess quick fix. this function does not work with changingVar either at the moment
	if (class(growObj)=="growthObjTruncIncr") { 
		predVar <- "incr"	
	} else {
		predVar <- attr(growObj@fit$terms,"predvars")[[2]]
	}
		
	if(predVar == "sizeNext") {
		plot(dataf$size, dataf$sizeNext,pch=19, xlab="Size at t", ylab="Size at t+1", main = mainTitle,...)
		abline(a = 0, b = 1)
	}else{
		dataf$incr <- dataf$sizeNext - dataf$size
		plot(dataf$size, dataf$incr, pch = 19, xlab = "Size at t", ylab="Size increment", main = mainTitle,...)
		abline(a = 0, b = 0)
	}
	
	if (length(grep("covariate", names(growObj@fit$model))) > 0) {  
		#convert to 1:n for indexing later and to relate to discrete
		dataf$covariate <- as.factor(dataf$covariate)
		levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
		ud <- unique(dataf$covariate); ud <- ud[!is.na(ud)]
		for (k in 1:length(ud)) { 
			tp <- dataf$covariate == ud[k]
			points(dataf$size[tp], dataf$sizeNext[tp], pch = 19, col = k)            
		}
		ud <- as.factor(ud)
	} else {
		ud <- 0
	} 
	
	sizes <- dataf$size[!is.na(dataf$size)]; sizes <- sizes[order(sizes)]
	
	for (k in 1:length(ud)) { 
		newd <- data.frame(size = sizes, size2 = sizes ^ 2,size3 = sizes ^ 3,
				covariate = as.factor(rep(ud[k],length(sizes))))
					
		if(length(grep("expsize", names(growObj@fit$coefficients))) == 1)
			newd$expsize = exp(sizes)
		if(length(grep("logsize", names(growObj@fit$coefficients))) == 1)
			newd$logsize = log(sizes)
		if(length(grep("logsize2", names(growObj@fit$coefficients))) == 1)
			newd$logsize=(log(sizes))^2
		
		
		if (length(grep("decline",tolower(as.character(class(growObj)))))>0 | 
				length(grep("trunc",tolower(as.character(class(growObj)))))>0) { 
				pred.size <- .predictMuX(growObj, newd, covPred = k)
			} else  {
				pred.size <- predict(growObj@fit,newd,type = "response")	
			}
		if (length(grep("incr", tolower(as.character(class(growObj))))) == 0) {
			points(sizes, pred.size, type = "l", col = k + 1)	
		} else { 
			if (length(grep("logincr", tolower(as.character(class(growObj))))) > 0) {
				points(sizes, sizes + exp(pred.size), type = "l", col = k + 1)		
			} else { 
				lines(sizes, pred.size, col = k + 1)	
			}
		}
	}
}


# =============================================================================
# =============================================================================
## FUNCTIONS FOR SIMULATING DATA ############################
## Return a data-frame with the headings used in the 'makeObject' functions


## Generate a simple data-frame for only continuous covariates
#  with a total of 1000 measurements with columns called
# size, sizeNext, surv, covariate, covariateNext, fec,
#
#
generateData <- function(nSamp=1000, type="simple"){
	if (nSamp<100) stop("Error: nSamp in generateData should be at least 100")
	if (type=="simple"){
	covariate <- sample(0:1, size=nSamp, replace=TRUE, prob = c(0.2, 0.8))
	covariateNext <- sample(0:1, size=nSamp, replace=TRUE, prob = c(0.8, 0.2))
	size <- rnorm(nSamp,5,2)
	#size <- exp(rnorm(1000, -1, 1.1))
	sizeNext <- 1+0.8*size-0.9*covariate+rnorm(nSamp,0,1)
	seedlings <- sample(1:nSamp,size=100,replace=TRUE)
	size[seedlings] <- NA; sizeNext[seedlings] <- rnorm(100,2,0.5)
	fec <- surv <- rep(NA, length(size))
	surv[!is.na(size)] <- rbinom(sum(!is.na(size)),1,invLogit(-1+0.2*size[!is.na(size)]))
	fec[!is.na(size)] <- rnorm(sum(!is.na(size)),exp(-7+0.9*size[!is.na(size)]),1)
	fec[size<quantile(size,0.20,na.rm=TRUE) | fec<0] <- 0
	fec <- fec*10
	
	stage <- stageNext <- rep("continuous",nSamp)
	stage[is.na(size)] <- NA
	stageNext[is.na(sizeNext)] <- "dead"
	
	dataf <- data.frame(size=size,sizeNext=sizeNext,surv=surv,
			covariate=covariate,covariateNext=covariateNext,
			fec=fec, stage=stage,stageNext=stageNext)
	
	dataf$sizeNext[dataf$surv==0] <- NA
	
	}
	
	if (type=="discrete") dataf <-.generateDataDiscrete(nSamp=nSamp)	
	if (type=="stochastic") dataf <-.generateDataStoch(nSamp=nSamp)	
	
	if (type!="simple" & type!="stochastic" & type!="discrete"){
		stop("unsupported type: supported types include simple, stochastic, or discrete")
	}
	
	return(dataf)
}

# =============================================================================
# =============================================================================
## FUNCTION FOR MAKING A LIST OF IPMS ############################################
# to do for stoch env with a single discrete covariate. ##########################

sampleSequentialIPMs <- function(dataf, nBigMatrix=10, minSize=-2,maxSize=10, 
		integrateType="midpoint", correction="none",
		explSurv=surv~size+size2+covariate,
		explGrow=sizeNext~size+size2+covariate, 
		regType="constantVar",explFec=fec~size,Family="gaussian", 
		Transform="none",fecConstants=data.frame(NA)) {
	
	#convert to 1:n for indexing later
	dataf$covariate <- as.factor(dataf$covariate)
	levels(dataf$covariate) <- 1:length(unique(dataf$covariate))
	
	print(explSurv)
	sv1 <- makeSurvObj(dataf=dataf,
			Formula=explSurv)
	gr1 <- makeGrowthObj(dataf=dataf,
			Formula=explGrow,
			regType=regType)
	
	fv1 <- makeFecObj(dataf=dataf,Formula=explFec, Family=Family, Transform=Transform, 
			fecConstants=fecConstants) 
	
	covs <- unique(dataf$covariate)
	covs <- covs[!is.na(covs)]
	
	#print(covs)
	
	IPM.list <- list()
	for (k in 1:length(covs)) { 
		
		tpF <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = data.frame(covariate=as.factor(k)),
				fecObj = fv1,integrateType=integrateType, correction=correction)
		tpS <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
				maxSize = maxSize, chosenCov = data.frame(covariate=as.factor(k)),
				growObj = gr1, survObj = sv1,
				integrateType=integrateType, correction=correction)
		IPM.list[[k]] <- tpF+tpS
	}
	return(IPM.list)
	
}


# =============================================================================
# =============================================================================
### simulation Carlina ######################################################

## Function to simulate something a bit like Carlina
##
## Parameters - nSamp - starting pop sizes
#             - nYrs - no years in simulation
#             - nSampleYrs - no years to extract for the data
#             - ... bunch of parameters
#             - meanYear - means for year effects
#             - matVarYear - variance covariances for year effects
#             - densDep - model density dependence in seedling establishment or not
#
# Returns - list including dataf - data-frame of data
#                                - various of the simulation parameters for convenience
#
#    PARAMETERS FROM TABLE 2 IN REES AND ELLNER 2009
simulateCarlina <- function(nSamp=200,nYrs=1000,nSampleYrs=15,
		m0=-1.37,ms=0.59,
		b0=-12.05,bs=3.64,
		A=-1,B=2,
		ag=1.14,bg=0.74,sig=0.29,
		mean.kids=3.16,sd.kids=0.5,
		meanYear=c(0,0,0),
		matVarYear=matrix(c(1.03,0,0,0,0.037,0.041,0,0.041,0.075),3,3),
		varA=0,varB=0,
		densDep=TRUE,maxPerYr=1000,maxStoreSeedlingsPerYr=200,sizes=c()) {
		
	require(mvtnorm)
	
	#initiate and set up year index
	if (length(sizes)==0) sizes <- rnorm(nSamp,3,0.5)
	startYr <- (nYrs-nSampleYrs)
	
	#recruits
	nrec <- c(20,42,12,17,8,19,58,45,44,2,56,25,75,92,94,6,4,34,104)
	
	#total plants
	totpl <- c(21,57,47,25,25,33,88,94,97,26,85,80,122,175,160,10,6,189)
		
	#set up dataframe
	dataf <-matrix(NA,(maxPerYr+maxStoreSeedlingsPerYr)*nYrs,13)
	maxPop <- nrow(dataf)
	n.per.yr  <-  rep(NA,nYrs)
	trueGrow <- rep(NA,nYrs)	
	
	count <- 0

	#stoch sims
	if (sum(matVarYear)!=0) {
		tmp <- rmvnorm(nYrs,mean=meanYear,sigma=matVarYear) 
	} else {
		tmp <- matrix(0,nYrs,3)
	}
	if (varA!=0) tmpA <- rnorm(nYrs,mean=A,sd=sqrt(varA))
	if (varB!=0) tmpB <- rnorm(nYrs,mean=B,sd=sqrt(varB))
	
	for (t in 1:nYrs) {
		
		#yr effects
		nSeedlings <- sample(nrec,size=1,replace=FALSE)
		#print(tmp)
		n.per.yr[t] <- length(sizes)

		m.year <- tmp[t,1]
		cg.year <- tmp[t,2]
		b.year <- tmp[t,3]
		
		if (varA!=0) A.year <- tmpA[t] else A.year <- A
		if (varB!=0) B.year <- tmpB[t] else B.year <- B
		
		
		#print(c(m.year,cg.year,b.year))
		
		if (n.per.yr[t]>0) { 
			#survival
			sx <- 1*(invLogit(m0+ms*sizes+m.year)>runif(n.per.yr[t]))
			#flowering
			fx <- 1*(invLogit(b0+bs*sizes)>runif(n.per.yr[t]))
			#fertility
			seedsx <- exp(A.year+B.year*sizes)*fx*sx
			#seedsx[seedsx>0] <- rpois(sum(seedsx>0),seedsx[seedsx>0])
		} else {
			sx <- fx <- seedsx <- c()
			print(c("extinct in year ", t))
		}
		
		if (densDep) pEst <- min(nSeedlings/max(sum(seedsx, na.rm=TRUE),1),1) else pEst <- 1
		
		babies <- rnorm(ceiling(pEst*max(sum(seedsx,na.rm=TRUE),1)),
				mean.kids+b.year,sd.kids) 
		#if (count>0) print(c("yr",t-startYr,pEst*max(sum(seedsx,na.rm=TRUE),1),length(babies)))
		#will end up with nrec babies at least in density dependent case
		if (length(babies)<nSeedlings & densDep) nSeedlings <- length(babies)
		
		#growth
		sizeNext <- rnorm(length(sizes),ag+bg*sizes+cg.year,sig)
		
		#remove dead or flowered
		fx[sx==0] <- NA
		sizeNext[sx==0 | fx==1] <- NA
		
		#storage
		if (t > startYr) {
			#print(t)
			#don't 'do it at all if too big (so as not to do only adults, etc)
			if ((count+length(sizes))>maxPop) { print("large pop size, breaking");break()}
			if ((count+length(sizes)+length(babies))>maxPop) { print("large pop size, breaking");break()}
			
			chs <- (count+1):(count+length(sizes))
			dataf[chs,1] <- sizes
			dataf[chs,2] <- sizeNext
			dataf[chs,3] <- sx
			dataf[chs,4] <- fx
			dataf[chs,5] <- seedsx
			dataf[chs,6] <- t
			dataf[chs,7] <- nSeedlings
			dataf[chs,8] <- m.year
			dataf[chs,9] <- cg.year
			dataf[chs,10] <- b.year
			dataf[chs,12] <- A.year
			dataf[chs,13] <- B.year
			
			count <- count + length(sizes)
						
			nbabes.store <- min(length(babies),maxStoreSeedlingsPerYr)
			chs <- (count+1):(count+nbabes.store)
			
			dataf[chs,2] <- babies[1:nbabes.store]
			dataf[chs,6] <- t
			dataf[chs,7] <- nSeedlings
			dataf[chs,8] <- m.year
			dataf[chs,9] <- cg.year
			dataf[chs,10] <- b.year
			dataf[chs,11] <- "sexual"
			dataf[chs,12] <- A.year
			dataf[chs,13] <- B.year
			
			count <- count + nbabes.store
					
			#print(nbabes.store)	
			
			} 
		
						
			#new pop - i.e. those that grew, survive, and didn't flower; as well as babies
			sizes <- c(sizeNext[sx==1 & fx==0 & !is.na(fx)],babies)
			if (length(sizes)==0) print("extinct")

			#get true growth rate 
			trueGrow[t] <- log(length(sizes)/n.per.yr[t])
			
			#cull size at t t down to maxPerYr, if not density dependent
			if (!densDep) sizes <- sample(sizes,size=maxPerYr, replace=TRUE)

	}
		

	list.par <- list(m0=m0,ms=ms,
			b0=b0,bs=bs,
			A=A,B=B,varB=varB,
			ag=ag,bg=bg,sig=sig,
			mean.kids=mean.kids,sd.kids=sd.kids,
			meanYear=meanYear,
			matVarYear=matVarYear,
			nrec=nrec)
	
	dataf <- dataf[1:count,]
	
	dataf <- data.frame(dataf,stringsAsFactors = FALSE)	
	#print(table(dataf[,6]))

	#plot(cumsum(trueGrow[startYr:nYrs])/(1:length(trueGrow[startYr:nYrs])),type="l")
	vartrueGrow <- var(trueGrow[startYr:nYrs],na.rm=TRUE) 
	meantrueGrow <- mean(exp(trueGrow[startYr:nYrs]),na.rm=TRUE)  
	trueGrow <- mean(trueGrow[startYr:nYrs],na.rm=TRUE)  
	#abline(h=trueGrow)
	
	colnames(dataf) <- c("size","sizeNext","surv","flower","fec",
			"year","nSeedlings","m.year","cg.year","b.year","offspringNext",
			"A.year","B.year")
	
	dataf$size <- as.numeric(dataf$size)
	dataf$sizeNext <- as.numeric(dataf$sizeNext)
	dataf$surv <- as.numeric(dataf$surv)
	dataf$flower <- as.numeric(dataf$flower)
	dataf$fec <- as.numeric(dataf$fec)
	dataf$nSeedlings <- as.numeric(dataf$nSeedlings)
	dataf$year <- as.numeric(dataf$year)
	dataf$cg.year <- as.numeric(dataf$cg.year)
	dataf$m.year <- as.numeric(dataf$m.year)
	dataf$b.year <- as.numeric(dataf$b.year)	
	dataf$A.year <- as.numeric(dataf$A.year)	
	dataf$B.year <- as.numeric(dataf$B.year)	
	
	#print(table(dataf$year))
	
	dataf$fec[dataf$fec==0] <- NA
	
	return(list(dataf=dataf,meanYear=meanYear,matVarYear=matVarYear,
					list.par=list.par,trueGrow=trueGrow,
					vartrueGrow=vartrueGrow,meantrueGrow=meantrueGrow))
	
}

# =============================================================================
# =============================================================================
## Convert Increment - where exact dates of census vary but some multiplier of yearly increments
## are desired; this function takes a data-frame (with columns size, sizeNext,
## and, importantly exactDate, exactDatel)
## and returns a data-frame with sizeNext modified proportional to the time
## elapsed in the desired yaerly increments, adding an additional column denoted 'increment'. 
#
# Parameters - dataf - a dataframe with headings size, sizeNext, exactDate,exactDatel
#            - nYrs - the number of years between censuses desired (e.g. for CTFS data, 5 years)
#
# Returns - a dataframe with the same headings
convertIncrement <- function(dataf, nYrs=1) {
	incr <- dataf$sizeNext - dataf$size
	if (class(dataf$exactDatel) == "Date")  {
		timeElapsed <- (difftime(dataf$exactDatel,dataf$exactDate)[[1]])/(365*nYrs)
	} else {
		timeElapsed <- (dataf$exactDatel-dataf$exactDate)/(365*nYrs)
	}
	incrNew <- incr / timeElapsed
	dataf$sizeNext <- dataf$size + incrNew
	dataf$incr <- dataf$incrNew
	return(dataf)
	
}

# =============================================================================
# =============================================================================
# Function to compare model fits for growth and survival objects built with different linear combinations of covariates. 
# Growth can have multiple response forms. Uses .makeCovDf.  Separate plot functions. 
#
#
# Returns - list with a summary table of covariates and scores, and a sub-list of growth or survival objects.
#
growthModelComp <- function(dataf, 
		expVars = c(sizeNext~1, sizeNext~size, sizeNext~size + size2), 
		regressionType = "constantVar",
		testType = "AIC",
		makePlot = FALSE,
		mainTitle = "",
		plotLegend = TRUE,
		legendPos = "topright", ...) {
	varN <- length(expVars)
	typeN <- length(regressionType)
	treatN <- varN * typeN
	summaryTable <- data.frame()
	grObj <- vector("list", length = treatN)
	i <- 1
	for(v in 1:varN) {
		for(t in 1:typeN) {
			grObj[[i]] <- makeGrowthObj(dataf = dataf, regType = regressionType[t], Formula = expVars[[v]]) 
			if (length(grep("decline",tolower(as.character(class(grObj[[i]])))))>0) { 
				summaryTable <- rbind(summaryTable, cbind(as.character(unlist(expVars))[v], regressionType[t],  match.fun(testType)(grObj[[i]]@fit$fit)))
			} else { 
			summaryTable <- rbind(summaryTable, cbind(as.character(unlist(expVars))[v], regressionType[t],  match.fun(testType)(grObj[[i]]@fit)))
			}	
			i <- i + 1
		}
	}
	summaryTable <- as.data.frame(summaryTable)
	names(summaryTable) <- c("Exp. Vars", "Reg. Type", testType)
	outputList <- list(summaryTable = summaryTable, growthObjects = grObj)
	
	# PLOT SECTION #
	if(makePlot == TRUE) {
		plotGrowthModelComp(grObj = grObj, summaryTable = summaryTable, dataf = dataf, expVars = expVars,  testType = "AIC",  plotLegend = plotLegend, mainTitle = mainTitle, legendPos, ...)
	}
	return(outputList)
}

# =============================================================================
# =============================================================================
survModelComp <- function(dataf, 
		expVars = c(surv~1, surv~size, surv~size + size2), 
		testType = "AIC",
		makePlot = FALSE,
		mainTitle = "", ncuts = 20,
		plotLegend = TRUE, 
		legendPos = "bottomleft", ...) {
	varN <- length(expVars)
	treatN <- varN
	summaryTable <- data.frame()
	svObj <- vector("list", length = treatN)
	i <- 1
	for(v in 1:varN) {
		svObj[[i]] <- makeSurvObj(dataf = dataf, Formula = expVars[[v]]) 
		summaryTable <- rbind(summaryTable, cbind(as.character(unlist(expVars))[v], match.fun(testType)(svObj[[i]]@fit)))
		i <- i + 1
	}
	summaryTable <- as.data.frame(summaryTable)
	names(summaryTable) <- c("Exp. Vars", testType)
	outputList <- list(summaryTable = summaryTable, survObjects = svObj)
	
	# PLOT SECTION #
	if(makePlot == TRUE) {
		## this is the surv picture    
		plotSurvModelComp(svObj = svObj, summaryTable = summaryTable, dataf = dataf, expVars = expVars, testType = "AIC", plotLegend = plotLegend, mainTitle = mainTitle, ncuts = ncuts, legendPos = legendPos, ...)
	}
	return(outputList)
}	

# =============================================================================
# =============================================================================
# Plot functions for model comparison.  Plots the series of fitted models for growth and survival objects.  
# Can plot a legend with the model covariates and model test criterion scores (defaults to AIC).

plotGrowthModelComp <- function(grObj, summaryTable, dataf, expVars,  testType = "AIC", plotLegend = TRUE, mainTitle = "", legendPos = "topright", ...) {
	treatN <- length(grObj)
	sizeSorted <- unique(sort(dataf$size))
	if(length(grep("sizeNext", unlist(as.character(expVars))))) {
		y.lab <- "Size at t + 1"
		dataSizeNext <- dataf$sizeNext
	}
	if(length(grep("incr", unlist(as.character(expVars))))) {
		y.lab <- "Growth"
		dataSizeNext <- dataf$sizeNext - dataf$size
	}
	if(length(grep("logincr", unlist(as.character(expVars))))){
		y.lab <- "log(growth)"
		dataSizeNext <- log(dataf$sizeNext - dataf$size)
	}
	plot(dataf$size, dataSizeNext, pch = 19, xlab = "Size at t", ylab = y.lab, main = mainTitle, cex = 0.8,...)
	for(p in 1:treatN) {
		
		#PROBLEM HERE 
		expVarsHere <- paste(attr(terms((expVars[[p]])),"term.labels"),collapse="+")
	
		newd <- .makeCovDf(sizeSorted, expVarsHere)
		if (length(grep("decline",tolower(as.character(class(grObj[[p]])))))>0) {
			pred.size <- .predictMuX(grObj[[p]], newd)
		} else { pred.size <- predict(grObj[[p]]@fit, newd, type = "response")}
		lines(sizeSorted, pred.size, type = "l", col = (p + 1))
	}
	if(plotLegend) {
		legend(legendPos, legend = sprintf("%s: %s = %.1f", expVars, testType, as.numeric(as.character(summaryTable[,3]))), col = c(2:(p + 1)), lty = 1, xjust = 1, bg = "white")
	}
}

# =============================================================================
# =============================================================================
plotSurvModelComp <- function(svObj, summaryTable, dataf,  expVars, testType = "AIC", plotLegend = TRUE, mainTitle = "", ncuts = 20, legendPos = "bottomleft", ...) {
	treatN <- length(svObj)
	#ncuts <- 20  # survival bins
	os <- order(dataf$size)  # order size
	osSurv <- (dataf$surv)[os] # order survival data according to size
	osSize<-(dataf$size)[os] # ordered size data
	binnedSize <- tapply(osSize, as.numeric(cut(osSize, breaks=ncuts)), mean, na.rm = TRUE); # bin Size data
	binnedSurv <- tapply(osSurv, as.numeric(cut(osSize, breaks=ncuts)), mean, na.rm = TRUE) #bin Survival probabilities
	plot(binnedSize, binnedSurv, pch = 19, xlab = "Size at t", ylab = "Survival to t + 1", main = mainTitle, cex = 0.8,...)
	for(p in 1:treatN) {
		expVarsHere <- paste(attr(terms(expVars[[p]]),"term.labels"),collapse="+")
		newd <- .makeCovDf(osSize, expVarsHere[p])
		lines(dataf$size[order(dataf$size)], surv(dataf$size[os], data.frame(covariate=1), svObj[[p]]), col = (p + 1))           
	}
	if(plotLegend) {
		legend(legendPos, legend = sprintf("%s: %s = %.1f", expVars, testType, as.numeric(as.character(summaryTable[,2]))), col = c(2:(p + 1)), lty = 1, xjust = 1, bg = "white")
	}
}

# =============================================================================
# =============================================================================
## Function to add pdf to model comparison pics
addPdfGrowthPic <- function(respType = "sizeNext",
		sizesPlotAt=c(20,50,60),
		sizeRange=c(20,400),
		incrRange=c(-10,50), 
		scalar=100,
		growthObjList,
		cols=1:5,
		cov=data.frame(covariate=1),
		minShow=1e-2,
		jitt=2,  #how far apart should pdfs be if you are plotting several
		...){
	
	nval <- length(growthObjList)
	sizes <- seq(sizeRange[1],sizeRange[2],length=500)
	incr <- seq(incrRange[1],incrRange[2],length=500)
	logincr <- log(incr); logincr[!is.finite(logincr)] <- NA
	
	for (j in 1:nval) {
		#print(j)
		for (k in 1:length(sizesPlotAt)) {           
			if (respType=="sizeNext") {
				pred <- growth(sizesPlotAt[k],sizes,cov,growthObjList[[j]])*scalar
				pred[pred<minShow] <- NA
				points(sizesPlotAt[k]+pred+jitt*(j-1),
						sizes,type="l", col=cols[j],...)                
			}
			if (respType=="incr") {
				pred <- growth(sizesPlotAt[k],sizesPlotAt[k]+incr,cov,growthObjList[[j]])*scalar
				pred[pred<minShow] <- NA
				points(sizesPlotAt[k]+pred+jitt*(j-1),
						incr,type="l", col=cols[j],...)
			}
			if (respType=="logincr") {
				pred <- growth(sizesPlotAt[k],sizesPlotAt[k]+logincr,cov,growthObjList[[j]])*scalar
				pred[pred<minShow] <- NA
				#print(range(pred,na.rm=TRUE))
				points(sizesPlotAt[k]+pred+jitt*(j-1),
						logincr,type="l", col=cols[j],...)
			}
		}}
}

# =============================================================================
# =============================================================================
## Function to coerce Growth object to parameters and variance desired
coerceGrowthObj <- function(growthObj, coeff, sd){

	if (length(growthObj@fit$coefficients) !=length(coeff)) {
		print("warning: number of desired coefficients to not match number of growth object coefficients")
	} 	
	growthObj@fit$coefficients[] <- as.numeric(coeff)
	growthObj@sd[] <- as.numeric(sd)
	return(growthObj)
}

# =============================================================================
# =============================================================================
## Function to coerce Survival object to parameters desired
coerceSurvObj <- function(survObj, coeff){
	
	if (length(survObj@fit$coefficients) !=length(coeff)) print("warning: number of desired coefficients to not match number of growth object coefficients")
	survObj@fit$coefficients[] <- as.numeric(coeff)
	
	return(survObj)
}

# =============================================================================
# =============================================================================
## Parametric boostrap of vital rate objects
# created by cory on 6-4-13 as a modification of getListRegObjs()

sampleVitalRateObj <- function(Obj,nSamp=100,nDiscreteGrowthTransitions=NULL,nDiscreteOffspringTransitions=NULL,nOffspring=NULL) {
  require(mvtnorm)
  require(MCMCpack)
  
  objList <- list()
  #generate new set parameters from mvn
  
  # what is the difference between offspring splitter transitions and discrete
  # for discrete transition matrices objects
  if( sum(class(Obj)==c("discreteTrans","discreteTransInteger")) ) {
    if(is.null(nDiscreteGrowthTransitions)){ 
      stop('Error: Please specify the total number of transitions that were used to estimate the discrete transitions in Obj@discreteTrans to allow the function to estimate the appropriate variance for resampling them')
    } 
    if(!is.null(nDiscreteGrowthTransitions)){ 
      for (j in 1:nSamp) {
        objList[[j]] <- Obj
        for (k in 1:ncol(Obj@discreteTrans)){
          objList[[j]]@discreteTrans[,k] <- rdirichlet(1, nDiscreteGrowthTransitions * as.numeric(Obj@discreteTrans[,k]))
        }
      }
    }
  }
  
  
  # for growth and survival objects
  if( !sum(class(Obj)==c("fecObj","fecObjInteger","discreteTrans", "discreteTransInteger")) ) {
    newpar <- rmvnorm(nSamp, mean = Obj@fit$coefficients, sigma = vcov(Obj@fit))
    for (j in 1:nSamp) {
      objList[[j]] <- Obj
      objList[[j]]@fit$coefficients <- newpar[j,]
    }
  }
  
  # for fecundity objects
  if( sum(class(Obj)==c("fecObj","fecObjInteger")) ) {
    for (j in 1:nSamp) {
      
      # sample fecundity regression parameters
      for (k in 1:length(Obj@fitFec)) { 
        newpar <- rmvnorm(1, mean = Obj@fitFec[[k]]$coefficients, 
            sigma = vcov(Obj@fitFec[[k]]))
        objList[[j]] <- Obj
        objList[[j]]@fitFec[[k]]$coefficients <- newpar
      }
      
      # sample discrete transitions from dirichlet
      if(ncol(Obj@offspringSplitter)>1){ # only if there are discrete fec stages
        if(is.null(nDiscreteOffspringTransitions)){ 
          stop('Error: Please specify the total number of discrete transitions that were used to estimate the fecundity transition in Obj@offspringSplitter to allow the function to estimate the appropriate variance for resampling them')
        } 
        if(!is.null(nDiscreteOffspringTransitions)){
          objList[[j]]@offspringSplitter[] <- rdirichlet(1, nDiscreteOffspringTransitions * as.numeric(Obj@offspringSplitter))
        }
      }
      
      # sample the offspring sizes
      # TODO: OffspringRel currently does not store all the information from the regression, which it should so that this works with offspring models that have more than an intercept
      if(is.null(nOffspring)){ 
        stop('Error: Please specify the total number of offspring that were used to estimate the offspring size distribtuion in Obj@offspringRel to allow the function to estimate the appropriate variance for resampling them')
      } 
      # TODO: Make work with the new OffspringObjs	
      if(!is.null(nOffspring)){ 
        newpar3 <- rnorm(1,coefficients(Obj@offspringRel)[1], Obj@sdOffspringSize/sqrt(nOffspring))
        objList[[j]]@offspringRel$coefficients[] <- newpar3
        newpar4 <- rnorm(1,Obj@sdOffspringSize, sqrt(2*Obj@sdOffspringSize^4)/(nOffspring-1))
        objList[[j]]@sdOffspringSize <- newpar4
      }
    }	
  } 
  
  return(objList)
}

# =============================================================================
# =============================================================================
## Use list of vital rate objects to make List of IPMS
# created by cory on 6-4-13 as a modification of getListIPMs()

# TODO: MAKE WORK FOR OFFSPRING OBJS

sampleIPM<- function(
    growObjList=NULL,survObjList=NULL,fecObjList=NULL,
    offspringObjList=NULL, discreteTransList=1, 
    nBigMatrix,minSize,maxSize, 
    covariates=FALSE,envMat=NULL,
    integrateType="midpoint",correction="none",warn=TRUE) {
  
  if(!is.null(offspringObjList)){
    stop('Sorry, lists of offspring objects are not yet supported. However, specify offspring size distributions via makeFecObj() is supported. Please set offspringObjList=NULL, or modify the function to accept offspring objects and email it to us.')
  }
  # make the same number of samples for each vital rate object
  n.samples <- max(c(length(growObjList),length(survObjList),length(fecObjList)))
  
  if(!is.null(growObjList) & length(growObjList)<n.samples){
    if(warn) warning('Length of growth object list is less than the length of another vital rate object list, so some members of the growth object list have been repeated.')
    growthObjList <- sample(growObjList,size=n.samples,replace=T)
  }
  
  if(!is.null(survObjList) & length(survObjList)<n.samples ){
    if(warn) warning('Length of survival object list is less than the length of another vital rate object list, so some members of the survival object list have been repeated.')
    survObjList <- sample(survObjList,size=n.samples,replace=T)
  }
  
  if(!is.null(fecObjList) & length(fecObjList)<n.samples ){
    if(warn) warning('Length of fec object list is less than the length of another vital rate object list, so some members of the fec object list have been repeated.')
    fecObjList <- sample(fecObjList,size=n.samples,replace=T)
  }
  
  if(!is.null(growObjList) & length(discreteTransList)<n.samples ){
    # if(warn) warning('Length of discreteTrans list is less than the length of another vital rate object list, so some members of the discreteTrans list have been repeated.')
    discreteTransList <- sample(discreteTransList,size=n.samples,replace=T)
  }
  
  if(!is.null(growObjList)){
    PmatrixList <- .makeListPmatrix(growObjList=growObjList, survObjList=survObjList, nBigMatrix=nBigMatrix, minSize=minSize, maxSize=maxSize, cov=covariates, envMat=envMat,discreteTransList=discreteTransList, integrateType=integrateType, correction=correction) 
    matrixList <- PmatrixList
  }
  
  if(!is.null(fecObjList)){
    FmatrixList <- .makeListFmatrix(fecObjList, nBigMatrix=nBigMatrix, minSize=minSize, maxSize=maxSize, cov=covariates, envMat=envMat, integrateType=integrateType, correction=correction) 
    matrixList <- FmatrixList
  }	
  
  if(!is.null(growObjList) & !is.null(fecObjList)){
    for(i in 1:n.samples){ 
      matrixList[[i]] <- PmatrixList[[i]] # to get all the Pmatrix slots
      matrixList[[i]]@.Data <- PmatrixList[[i]]+FmatrixList[[i]] 
    }
  }
  
  return(matrixList)
}



# =============================================================================
# =============================================================================
## Use list of IPMs to make List of IPM output
# created by cory on 6-4-13 as a modification of getIPMOutput() and getIPMOutputDirect()

sampleIPMOutput <- function(IPMList=NULL,PMatrixList=NULL,passageTimeTargetSize=c(), sizeToAgeStartSize=c(),sizeToAgeTargetSize=c(),warn=TRUE){
  # removed from old version of getIPMOutputDirect
  # # 		cov=FALSE, envMat=NULL,
  # # 		chosenCov = data.frame(covariate = 1),
  # # 		onlyLowerTriGrowth=FALSE)
  
  # if only a Pmatrix is supplied, just call it IPMList for simplicity
  if(is.null(IPMList))  IPMList=PMatrixList
  if ( is.null(passageTimeTargetSize) )  { 
    if(warn) print("no target size for passage time provided; taking meshpoint median")
    passageTimeTargetSize <- median(IPMList[[1]]@meshpoints)
  }
  if ( is.null(sizeToAgeStartSize) )  { 
    if(warn) print("no starting size for size to age provided; taking minimum size")
    sizeToAgeStartSize <- min(IPMList[[1]]@meshpoints)
  }
  if ( is.null(sizeToAgeTargetSize) )  { 
    if(warn) print("no target size for size to age provided; taking meshpoint values")
    sizeToAgeTargetSize <- IPMList[[1]]@meshpoints
  }
  nSamps <- length(IPMList)
  
  stableStage <- LE <- pTime <- matrix(NA,nSamps,length(IPMList[[1]]@.Data[1,]))
  lambda <- rep(NA,nSamps)
  resAge <- resSize <- matrix(NA,nSamps,length(sizeToAgeTargetSize)) 
  h1 <- diff(IPMList[[1]]@meshpoints)[1]
  
  for (k in 1:nSamps) {
    IPM <- IPMList[[k]]
    LE[k,]<-meanLifeExpect(IPM) 
    pTime[k,]<-passageTime(passageTimeTargetSize,IPM) 
    
    if (is.null(PMatrixList)) {
      lambda[k] <- Re(eigen(IPM)$value[1])
      stableStage[k,] <- eigen(IPM)$vector[,1]
      #normalize stable size distribution
      stableStage[k,] <- stableStage[k,]/sum(stableStage[k,])
    }
    # get size to age results
    res2 <- sizeToAge(Pmatrix=IPM,startingSize=sizeToAgeStartSize, 
        targetSize=sizeToAgeTargetSize)
    resAge[k,] <- res2$timeInYears
    resSize[k,] <- res2$targetSize
  }
  
  return(list(LE=LE,pTime=pTime,lambda=lambda,stableStage=stableStage,
          resAge=resAge,resSize=resSize,meshpoints=IPM@meshpoints))	
}

# =============================================================================
# =============================================================================
invLogit <- function(x) { 
  u <- exp(pmin(x, 50))
  return(u / (1 + u))
}

# =============================================================================
# =============================================================================
# For a single Pmatrix (!not compound and no discrete stages!), this functions makes a series
# of diagnostic plots - this is defined for growthObj,
# growthObjIncr - modification required
# if other objects used
#
# Parameters - the Pmatrix
#            - growObj - growth object used to build it
#            - survObj - survival object used to build it
#            - dff - the data from which it was built
#            - integrateType - "midpoint", or "cumul" - should
#                 correspond to what the IPM was built with
#            - do you want to implement the corrections for integration? 
# Returns - 
#

diagnosticsPmatrix <- function (Pmatrix, growObj, survObj, dff=NULL, 
    integrateType = "midpoint", 
    correction = "none", cov = data.frame(covariate = 1), 
    sizesToPlot=c(),extendSizeRange=c()) {
  print("Range of Pmatrix is ")
  print(range(c(Pmatrix)))
  if (Pmatrix@meshpoints[1] > 0) 
    new.min <- Pmatrix@meshpoints[1]/2
  else new.min <- Pmatrix@meshpoints[1] * 1.5
  new.max <- 1.5 * max(Pmatrix@meshpoints)
  
  if (length(extendSizeRange)>0 & length(extendSizeRange)!=2) print("require two values for extendSizeRange, reflecting upper and lower limits")
  if (length(extendSizeRange)>0) { new.min <- extendSizeRange[1]; new.max <- extendSizeRange[2]}
  
  
  #colours for 1) current; 2) bigger size range; 3) bigger no bins
  cols <- c("black","tomato","darkblue")
  ltys <- c(1,1,3)
  
  #matrix with bigger size range
  Pmatrix1 <- makeIPMPmatrix(nEnvClass = 1, nBigMatrix = length(Pmatrix@meshpoints), 
      minSize = new.min, maxSize = new.max, 
      chosenCov = cov, growObj = growObj, survObj = survObj, 
      integrateType = integrateType, correction = correction)
  
  if (sum(is.na(Pmatrix1)) > 0) {
    print("Pmatrix with extended size range returns NAs; changing these to 0, and putting survival value onto diagonal for columns that sum to zero")
    Pmatrix1[is.na(Pmatrix1)] <- 0
    bad <- which(colSums(Pmatrix1) == 0, arr.ind = TRUE)
    if (length(bad) > 0) 
      Pmatrix1[cbind(bad, bad)] <- surv(size = Pmatrix1@meshpoints[bad], 
          cov = cov, survObj = survObj)
  }
  
  
  #matrix with bigger number of bins	  
  Pmatrix2 <- makeIPMPmatrix(nEnvClass = 1, nBigMatrix = floor(length(Pmatrix@meshpoints) * 
              1.5), minSize =Pmatrix@meshpoints[1], maxSize = max(Pmatrix@meshpoints), 
      chosenCov = cov, growObj = growObj, survObj = survObj, 
      integrateType = integrateType, correction = correction)
  
  if (sum(is.na(Pmatrix2)) > 0) {
    print("Pmatrix with extended number of bins returns NAs; changing these to 0, and putting survival value onto diagonal for columns that sum to zero")
    Pmatrix2[is.na(Pmatrix2)] <- 0
    bad <- which(colSums(Pmatrix2) == 0, arr.ind = TRUE)
    if (length(bad) > 0) 
      Pmatrix2[cbind(bad, bad)] <- surv(size = Pmatrix2@meshpoints[bad], 
          cov = cov, survObj = survObj)
  }
  
  #start plots - put original Pmatrix in black
  #par(mfrow = c(3, 3), bty = "l")
  
  #save plot characters
  old.par <- par(no.readonly = TRUE) 
  
  par(mfrow = c(1, 3), bty = "l",pty="s", mar=c(5,4,4,1))
  if (!is.null(dff)) { 
    xlims <- range(c(Pmatrix@meshpoints, Pmatrix1@meshpoints,
            dff$size, dff$sizeNext), na.rm = TRUE)
    a1 <- hist(c(dff$size, dff$sizeNext), xlab = "Sizes observed", axes=FALSE,
        ylab = "Frequency", main = "", xlim = xlims, col="lightgrey", border="lightgrey",plot=TRUE)
    axis(1); axis(2)
  } else {
    a1 <- list(); a1$counts <- 1:100
    xlims <- range(c(Pmatrix@meshpoints, Pmatrix1@meshpoints))
    plot(1:100,type="n",xlab = "Sizes", ylab="", ylim=range(a1$counts), xlim=xlims)	
  } 
  
  lcs <- c(0.8,0.6,0.4)	
  lcs.x <- mean(xlims) 
  
  text(lcs.x ,max(a1$counts)*lcs[1],"Current",pos=3,col=cols[1])
  arrows(Pmatrix@meshpoints[1], max(a1$counts)*lcs[1],Pmatrix@meshpoints[length(Pmatrix@meshpoints)], max(a1$counts)*lcs[1],col=cols[1], length=0.1, code=3,lty=ltys[1])
  text(lcs.x ,max(a1$counts)*lcs[2],"Extended range",pos=3,col=cols[2])
  arrows(Pmatrix1@meshpoints[1], max(a1$counts)*lcs[2],Pmatrix1@meshpoints[length(Pmatrix1@meshpoints)], max(a1$counts)*lcs[2],col=cols[2], length=0.1, code=3,lty=ltys[1])
  text(lcs.x ,max(a1$counts)*lcs[3],"Increased no of bins",pos=3,col=cols[3])
  arrows(Pmatrix2@meshpoints[1], max(a1$counts)*lcs[3],Pmatrix2@meshpoints[length(Pmatrix2@meshpoints)], max(a1$counts)*lcs[3],col=cols[3], length=0.1, code=3,lty=ltys[1])
  
  #legend("topright", legend = "fitted range", col = "black", lty = 2, bty = "n") ##!!! change this
  title("Size range")
  
  #survival sums
  lims <- range(c(colSums(Pmatrix), surv(Pmatrix@meshpoints, cov, survObj)))
  plot(colSums(Pmatrix), surv(Pmatrix@meshpoints, cov, survObj), 
      type = "n", xlab = "Surviving in Pmatrix", ylab = "Should be surviving",col=cols[1],lty=ltys[1], 
      xlim=lims,ylim=lims)
  abline(0, 1, lwd=2,col="grey")
  points(colSums(Pmatrix), surv(Pmatrix@meshpoints, cov, survObj), type = "l", col = cols[1],lty=ltys[1])
  points(colSums(Pmatrix1), surv(Pmatrix1@meshpoints, cov, survObj), type = "l", col = cols[2],lty=ltys[2])
  points(colSums(Pmatrix2), surv(Pmatrix2@meshpoints, cov, survObj), type = "l", col = cols[3],lty=ltys[3])
  title("Survival")
  #text(lims[2],lims[2],"(0,1)",pos=1)	
  
  # mtext("Points should lie along the 0,1 line, shown in grey", side=1,outer=TRUE,line=-1)
  
  
  
  LE <- meanLifeExpect(Pmatrix)
  LE1 <- meanLifeExpect(Pmatrix1)
  LE2 <- meanLifeExpect(Pmatrix2)
  
  plot(Pmatrix@meshpoints, LE, type = "l", xlim = range(Pmatrix1@meshpoints), 
      ylim = range(c(LE, LE1)), xlab = "Sizes", ylab = "Life expectancy",col=cols[1],lty=ltys[1])
  points(Pmatrix1@meshpoints, LE1, type = "l", col = cols[2],lty=ltys[2])
  points(Pmatrix2@meshpoints, LE2, type = "l", col = cols[3],lty=ltys[3])
  legend("topleft", legend = c("Current", "Extended range", "Increased no of bins"), col = cols, lty = ltys, bty = "n")
  title("Life expectancy")
  
  #mtext("Increasing size range (red) or number of bins (blue) should not alter life expectancy estimates", side=1,outer=TRUE,line=-1)
  
  
  
  if (length(sizesToPlot)==0 & !is.null(dff)) sizesToPlot <- as.numeric(quantile(dff$size,c(0.25, 0.5,0.75),na.rm=TRUE))
  if (length(sizesToPlot)==0 & is.null(dff)) sizesToPlot <- as.numeric(quantile(Pmatrix@meshpoints,c(0.25, 0.5,0.75),na.rm=TRUE))
  
  loctest <- rep(NA,length(sizesToPlot))
  
  print("Please hit any key for the next plot")	
  scan(quiet="TRUE")
  
  par(mfrow = c(2, 3), bty = "l",pty="s")
  
  
  for (kk in c(1,3)) {
    if (kk==1) Pmat <- Pmatrix	
    if (kk==3) Pmat <- Pmatrix2	
    
    h <- diff(Pmat@meshpoints)[1]
    testSizes <- seq(min(Pmat@meshpoints), max(Pmat@meshpoints), 
        length = 5000)
    for (j in 1:3) {
      
      loctest[j] <- which(abs(Pmat@meshpoints-sizesToPlot[j])==min(abs(Pmat@meshpoints-sizesToPlot[j])),arr.ind=TRUE)[1]
      
      #print(loctest[j])
      
      ps <- surv(Pmat@meshpoints[loctest[j]], cov, survObj)
      newd <- data.frame(size = Pmat@meshpoints[loctest[j]], 
          size2 = Pmat@meshpoints[loctest[j]]^2, 
          size3 = Pmat@meshpoints[loctest[j]]^3, 
          covariate = Pmat@env.index[1])
      if (length(grep("expsize", names(growObj@fit$coefficients)))) 
        newd$expsize = log(Pmat@meshpoints[loctest[j]])
      if (length(grep("logsize", names(growObj@fit$coefficients)))) 
        newd$logsize = log(Pmat@meshpoints[loctest[j]])
      if (length(growObj@fit$model$covariate) > 0) 
        if (is.factor(growObj@fit$model$covariate)) 
          newd$covariate <- as.factor(newd$covariate)
      if (length(grep("decline", tolower(as.character(class(growObj))))) > 0 | 
          length(grep("trunc", tolower(as.character(class(growObj))))) > 0) {
        mux <- .predictMuX(growObj, newd)
      } else {
        mux <- predict(growObj@fit, newd, type = "response")
      }
      if (length(grep("incr", tolower(as.character(class(growObj))))) >  0 & 
          length(grep("logincr", tolower(as.character(class(growObj))))) == 0) {
        mux <- Pmat@meshpoints[loctest[j]] + mux
      }
      if (length(grep("decline", tolower(as.character(class(growObj))))) ==  0 &
          length(grep("trunc", tolower(as.character(class(growObj))))) == 0) {
        sigmax2 <- growObj@sd^2
      } else {
        sigmax2 <- growObj@fit$sigmax2
        var.exp.coef <- growObj@fit$var.exp.coef
        sigmax2 <- sigmax2 * exp(2 * (var.exp.coef * mux))
        if (length(grep("trunc", tolower(as.character(class(growObj))))) > 0) 
          sigmax2 <- growObj@fit$sigmax2
      }
      
      #  range.x <- range(Pmat@meshpoints[loctest[j]] + c(-3.5 * sqrt(sigmax2), +3.5 * sqrt(sigmax2)))
      range.x <- range(mux + c(-3.5 * sqrt(sigmax2), +3.5 * sqrt(sigmax2)))
      
      plot(Pmat@meshpoints, Pmat@.Data[, loctest[j]]/h/ps, 
          type = "n", xlim = range.x, xlab = "Size next", ylab = "Kernel")
      
      #if (j == 1 & kk==1) title("Numerical resolution and growth")
      if (j == 2 & kk==1) mtext("Numerical resolution and growth",3,line=4,font=2)
      if (j ==1 & kk==1) title("Current Pmatrix")
      if (j ==1 & kk==3) title("Increased bins")
      
      
      
      for (k in 1:length(Pmat@meshpoints)) {
        points(c(Pmat@meshpoints[k]) + c(-h/2, h/2), rep(Pmat@.Data[k,loctest[j]], 2)/h/ps, type = "l",col=cols[kk])
        points(rep(Pmat@meshpoints[k] + h/2, 2), c(0, Pmat@.Data[k, loctest[j]]/h/ps), type = "l", 
            lty = 1,col=cols[kk])
        points(rep(Pmat@meshpoints[k] - h/2, 2), c(0,Pmat@.Data[k, loctest[j]]/h/ps), type = "l", 
            lty = 1,col=cols[kk])
      }
      if (length(grep("logincr", tolower(as.character(class(growObj))))) == 0 &
          length(grep("trunc", tolower(as.character(class(growObj))))) == 0) {
        points(testSizes, dnorm(testSizes, mux, sqrt(sigmax2)), type = "l", col = 2)
      } else {
        if (length(grep("trunc", tolower(as.character(class(growObj))))) > 0) {
          require(truncnorm)
          points(testSizes, dtruncnorm(testSizes, a = Pmat@meshpoints[loctest[j]], 
                  b = Inf, mean = mux, sd = sqrt(sigmax2)), type = "l",col = 2)
        } else {
          points(testSizes, dlnorm(testSizes - Pmat@meshpoints[loctest[j]], 
                  mux, sqrt(sigmax2)), type = "l", col = 2)
        }
      }
      
      legend("topright",legend=paste("size=",round(Pmat@meshpoints[loctest[j]],2)), col = "white", 
          lty = 1, bty = "n")
      
    }
  }
  
  #reset plot characters
  on.exit(par(old.par))
  
}

# =============================================================================
# =============================================================================
# show the news file

IPMpackNews <- function () 
{
  newsfile <- file.path(system.file(package = "IPMpack"), 
      "IPMpackNews")
  file.show(newsfile)
}



