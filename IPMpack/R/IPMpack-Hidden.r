# TODO: Add comment
# 
# Author: ctg
###############################################################################
# This file contains functions that are hidden because they are no longer supported or because we like hiding them.

# =============================================================================
# =============================================================================
## Generate a simple data-frame for continuous and discrete covariates
#  with a total of 1000 measurements in columns called
# size, sizeNext, surv, fec, stage, stageNext number
# Stage contains names including "continuous", and then a range
# of names for discrete stages, e.g., in this example,
#  "dormant" "seedAge1"   "seedOld" 
#
# 
.generateDataDiscrete <- function(nSamp=1000){
  size <- rnorm(nSamp,5,2)
  sizeNext <- 1+0.8*size+rnorm(nSamp,0,1)
  surv <- rbinom(nSamp,1,invLogit(-1+0.2*size))
  sizeNext[surv==0] <- NA
  fec <- rnorm(length(size),exp(-7+0.9*size),1)
  fec[size<quantile(size,0.20) | fec<0] <- 0
  stage <- rep("continuous",nSamp)
  stageNext <- rep("continuous",nSamp)
  sizeNext[surv==0] <- NA
  stageNext[surv==0] <- c("dead")
  number <- rep(1,nSamp)
  become.dormant <- which(rank(size)%in%sample(rank(size),50,prob=surv*fec))
  sizeNext[become.dormant] <- NA; stageNext[become.dormant] <- c("dormant")
  were.dormant <- which(rank(sizeNext)%in%sample(rank(sizeNext),50,prob=surv*fec))
  size[were.dormant] <- NA; stage[were.dormant] <- c("dormant")
  dataf <- rbind(data.frame(size=size,sizeNext=sizeNext,surv=surv,
          fec=fec,stage=stage,stageNext=stageNext,number=number),
      data.frame(size=NA,sizeNext=NA,surv=rep(c(1,0),2),fec=0,
          stage=rep(c("seedAge1","seedOld"),each=2),
          stageNext=rep(c("seedOld","dead"),2),
          number=c(202,220,115,121)),
      data.frame(size=NA,sizeNext=rnorm(113,3,2),surv=1,fec=0,
          stage=c(rep("seedAge1",33),rep("seedOld",30),rep(NA,50)),
          stageNext=c("continuous"),number=1))
  
  
  return(dataf)
}

# =============================================================================
# =============================================================================
## Generate a simple data-frame for continuous and discrete covariates
#  with a total of 1000 measurements in columns called
# size, sizeNext, surv, fec, stage, stageNext number
#
# 
.generateDataStoch <- function(nSamp=1000){
  covariate1 <- rnorm(nSamp)
  covariate2 <- rnorm(nSamp)
  covariate3 <- rnorm(nSamp)
  size <- rnorm(nSamp,5,2)
  sizeNext <- 1+0.9*size+3*covariate1+0.01*covariate2+0.2*covariate3+rnorm(nSamp,0,0.1)
  
  fec <- surv <- rep(NA, length(size))
  surv[!is.na(size)] <- rbinom(sum(!is.na(size)),1,invLogit(-1+0.2*size[!is.na(size)]))
  fec[!is.na(size)] <- rnorm(sum(!is.na(size)),exp(-7+0.9*size[!is.na(size)]),1)
  fec[size<quantile(size,0.20,na.rm=TRUE) | fec<0] <- 0
  fec <- 10*fec
  
  seedlings <- sample(1:nSamp,size=100,replace=TRUE)
  size[seedlings] <- NA; 
  sizeNext[seedlings] <- rnorm(100,-2,0.1)
  surv[seedlings] <- 1
  #set to flower when covariate1 is around 1.5
  pfec <- 1*(runif(length(size))<invLogit(size+covariate1)); #print(pfec)
  fec[pfec==0] <- 0
  #fill in stage
  stage <- stageNext <- rep("continuous",nSamp)
  stage[is.na(size)] <- NA
  stageNext[is.na(sizeNext)] <- "dead"
  
  dataf <- data.frame(size=size,sizeNext=sizeNext,surv=surv,
      covariate1=covariate1,covariate2=covariate2,covariate3=covariate3,
      fec=fec, stage=stage,stageNext=stageNext, number=rep(1,length(size)))
  
  dataf$sizeNext[dataf$surv==0] <- NA
  
  return(dataf)
}

# =============================================================================
# =============================================================================
# Function to extract IPM output from a list
# of P (survival + growth) and F (fecundity) matrices
# (usually from Bayes fit) 
#
# Parameters - PmatrixList
#            - targetSize - the size you want passage time estimated for.
#            - FmatrixList
#
# Returns - a list 

.getIPMoutput <- function(PmatrixList,targetSize=c(),FmatrixList=NULL){
  
  if (length(targetSize)==0)  { 
    print("no target size for passage time provided; taking meshpoint median")
    targetSize <- median(PmatrixList[[1]]@meshpoints)
  }
  nsamps <- length(PmatrixList)
  h1 <- PmatrixList[[1]]@meshpoints[2]-PmatrixList[[1]]@meshpoints[1]
  stableStage <- LE <- pTime <- matrix(NA,nsamps,length(PmatrixList[[1]]@.Data[1,]))
  lambda <- rep(NA,nsamps)
  for (k in 1:nsamps) {
    Pmatrix <- PmatrixList[[k]]
    LE[k,]<-meanLifeExpect(Pmatrix) 
    pTime[k,]<-passageTime(targetSize,Pmatrix) 
    
    if (class(FmatrixList)!="NULL") {
      IPM <- Pmatrix + FmatrixList[[k]]
      lambda[k] <- Re(eigen(IPM)$value[1])
      stableStage[k,] <- eigen(IPM)$vector[,1]
      #normalize stable size distribution
      stableStage[k,] <- stableStage[k,]/(h1*sum(stableStage[k,]))
    }
  }
  
  return(list(LE=LE,pTime=pTime,lambda=lambda,stableStage=stableStage))
  
}

# =============================================================================
# =============================================================================
# Function to extract IPM output from posteriors 
# (usually from Bayes fit)  - quicker to do all at once
# rather than build list of IPM P matrices, then list of IPM F matrices
#
# Parameters - survObjlist - list of survival objects
#            - growObjList - list of growth objects
#            - targetSize - the size you want passage time estimated for.
#            - nBigMatrix - the number of bins
#            - minSize - the minimum size
#            - maxSize - the maximum size
#            - cov - do you want to fit a discrete covariate
#            - fecObjList - list of fecundity objects
#            - envMat - matrix of env transitions (only if cov=TRUE)
#            - nsizeToAge - numeric describing how many size to age defined (0 - 100s)
# 
#
# Returns - a list 

.getIPMoutputDirect <- function(survObjList,growObjList,targetSize=c(),
    nBigMatrix,minSize,maxSize,discreteTrans = 1,
    cov=FALSE,fecObjList=NULL, envMat=NULL,
    nsizeToAge=0, sizeStart = 10,
    integrateType = "midpoint", correction = "none", storePar=TRUE,
    chosenCov = data.frame(covariate = 1),
    onlyLowerTriGrowth=FALSE){
  
  # adjust the sample lengths to they are all the same
  if (length(targetSize)==0)  targetSize <- 0.2*(minSize+maxSize)
  nsamp <- max(length(growObjList),length(survObjList),length(fecObjList))
  if (length(survObjList)<nsamp)  
    survObjList <- sample(survObjList,size=nsamp,replace=TRUE)
  if (length(growObjList)<nsamp)  
    growObjList <- sample(growObjList,size=nsamp,replace=TRUE)
  if (class(fecObjList)!="NULL") {
    if (length(fecObjList)<nsamp)  
      fecObjList <- sample(fecObjList,size=nsamp,replace=TRUE)
  }
  
  # store chosen parameters
  if (storePar){
    surv.par <- matrix(NA,nsamp,length(survObjList[[1]]@fit$coefficients))
    grow.par <- matrix(NA,nsamp,length(growObjList[[1]]@fit$coefficients)+1)
    for (k in 1:nsamp) {
      surv.par[k,] <- survObjList[[k]]@fit$coefficients
      grow.par[k,] <- c(growObjList[[k]]@fit$coefficients,growObjList[[k]]@sd)
    }} else { surv.par <- grow.par <- c()}
  
  #set up storage
  if (class(discreteTrans)=="discreteTrans") nDisc <- (ncol(discreteTrans@discreteTrans)-1) else nDisc <- 0
  
  if (class(envMat)!="NULL") nEnv <- envMat@nEnvClass else nEnv <- 1
  LE <- pTime <- matrix(NA,nsamp,(nBigMatrix+nDisc)*nEnv)
  if (class(fecObjList)=="NULL") {
    lambda <- stableStage <- c()
  } else {
    stableStage <- matrix(NA,nsamp,(nBigMatrix+nDisc)*nEnv)
    lambda <- rep(NA,nsamp)
  }
  if (nsizeToAge==0) { resAge <- resSize <- c() } else { resAge <- resSize <- matrix(NA,nsamp,nsizeToAge)} 
  if (length(sizeStart)==0) { if (minSize<0) sizeStart <- 0.5*minSize else sizeStart <- 2*minSize }
  
  #go!
  for (k in 1:nsamp) {
    
    if (!cov) {
      Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
          maxSize = maxSize,  growObj = growObjList[[k]],
          survObj = survObjList[[k]],discreteTrans=discreteTrans,
          integrateType=integrateType, correction=correction) 
      
    } else {
      Pmatrix <- makeCompoundPmatrix(nEnvClass = nEnv,
          nBigMatrix = nBigMatrix, minSize = minSize,
          maxSize = maxSize, envMatrix=envMat,growObj = growObjList[[k]],
          survObj = survObjList[[k]],discreteTrans=discreteTrans,
          integrateType=integrateType, correction=correction)    
      
    }
    
    if (onlyLowerTriGrowth & !cov) {
      Pmatrix@.Data <- Pmatrix@.Data*lower.tri(Pmatrix@.Data, diag = TRUE)
      nvals <- colSums(Pmatrix@.Data,na.rm=TRUE)
      Pmatrix@.Data <- t((t(Pmatrix@.Data)/nvals) *
              surv(size = Pmatrix@meshpoints, 
                  cov = chosenCov,
                  survObj = survObjList[[k]]))                                    
    }
    
    
    LE[k,] <- meanLifeExpect(Pmatrix) 
    pTime[k,] <- passageTime(targetSize,Pmatrix) 
    if (k==1) h1 <- diff(Pmatrix@meshpoints)[1]
    
    if (class(fecObjList)!="NULL") {
      if (!cov) { 
        Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
            maxSize = maxSize,  
            fecObj=fecObjList[[k]],
            integrateType=integrateType, correction=correction)
      } else {
        Fmatrix <- makeCompoundFmatrix(nEnvClass = nEnv,
            nBigMatrix = nBigMatrix, minSize = minSize, 
            maxSize = maxSize, envMatrix=envMat,
            fecObj=fecObjList[[k]],integrateType=integrateType, correction=correction)
      }
      
      
      
      IPM <- Pmatrix + Fmatrix
      lambda[k] <- Re(eigen(IPM)$value[1])
      stableStage[k,] <- eigen(IPM)$vector[,1]
      #normalize stable size distribution
      stableStage[k,] <- stableStage[k,]/(h1*sum(stableStage[k,]))
      
      #print("here2")
    }
    
    # get size to age results
    if (nsizeToAge>0) { 
      res2 <- sizeToAge(Pmatrix=Pmatrix,startingSize=minSize*1.3,
          targetSize=seq(sizeStart,maxSize*0.9,length=nsizeToAge))
      resAge[k,] <- res2$timeInYears
      resSize[k,] <- res2$targetSize
    }
    
  }
  
  return(list(LE=LE,pTime=pTime,lambda=lambda,stableStage=stableStage,
          meshpoints=Pmatrix@meshpoints,resAge=resAge,resSize=resSize,
          surv.par=surv.par,grow.par=grow.par))
  
}

# =============================================================================
# =============================================================================
## Function to take fit of these and output a list of growth objects
.getListRegObjects <- function(Obj,nsamp=1000) {
  
  require(mvtnorm)
  
  #generate new set parameters from mvn
  npar <- length(Obj@fit$coefficients)
  newpar <- rmvnorm(nsamp, mean = Obj@fit$coefficients, sigma = vcov(Obj@fit))
  
  objList <- list()
  
  for (j in 1:nsamp) {
    objList[[j]] <- Obj
    objList[[j]]@fit$coefficients <- newpar[j,]
  }
  
  return(objList)
}


# =============================================================================
# =============================================================================
## Function to take fit of these and output a list of growth objects
.getListRegObjectsFec <- function(Obj,nsamp=1000) {
  
  require(mvtnorm)
  
  objList <- list()
  
  #generate new set parameters from mvn
  for (j in 1:nsamp) {
    for (k in 1:length(Obj@fitFec)) { 
      npar <- length(Obj@fitFec[[k]]$coefficients)
      newpar <- rmvnorm(nsamp, mean = Obj@fitFec[[k]]$coefficients, 
          sigma = vcov(Obj@fitFec[[k]]))
      objList[[j]] <- Obj
      objList[[j]]@fitFec[[k]]$coefficients <- newpar[j,]
    }
  }	
  
  return(objList)
}

# =============================================================================
# =============================================================================

#Find years where can estimate all three stochastic vital rates(survival, growth and baby size)
.identifyPossibleYearsCarlina <- function(dataf){
  
  yr1 <- table(dataf$year[!is.na(dataf$size) & 
              !is.na(dataf$sizeNext) & is.na(dataf$offspringNext)])
  yr2 <- table(dataf$year[!is.na(dataf$size) & 
              !is.na(dataf$surv) & is.na(dataf$offspringNext)])
  yr3 <- table(dataf$year[!is.na(dataf$sizeNext) & 
              !is.na(dataf$offspringNext)])
  
  good.yrs <- intersect(as.numeric(as.character(names(yr1)[yr1>1])),
      as.numeric(as.character(names(yr2))[yr2>1]))
  good.yrs <- intersect(good.yrs,as.numeric(as.character(names(yr3)[yr3>1])))
  
  return(is.element(dataf$year,good.yrs))
}

# =============================================================================
# =============================================================================
# Function to plot the results of a stochastic simulation
# structure run
#
# Parameters - tVals - time points
#            - st - output of stochGrowthRateManyCov
#            - covTest - the key covariate for germination / flowering
#            - nRunIn - how many to leave off pics
# Returns - 
#
.plotResultsStochStruct <- function(tVals,meshpoints,st,covTest,nRunIn=15,log="y",...) { 
  
  par(mfrow=c(1,3),bty="l", pty="s")
  plot(tVals[nRunIn:length(tVals)],
      colSums(st[,nRunIn:length(tVals)]),
      xlab="Time", 
      ylab="Population size",type="l",...)
  abline(v=1:max(tVals),lty=3)
  
  for (j in 1:ncol(covTest)) {
    #normalized
    covTest[,j] <- (covTest[,j]-mean(covTest[,j]))/sd(covTest[,j])
  }
  
  matplot(tVals[nRunIn:length(tVals)],covTest[nRunIn:length(tVals),],
      type="l",lty=3,col=1:ncol(covTest),xlab="Time", ylab="Covariates")
  abline(v=1:max(tVals),lty=3)
  
  
  if (log=="y") st <- log(st)
  
  image(tVals[nRunIn:length(tVals)],
      meshpoints,
      t(st[,nRunIn:length(tVals)]),
      ylab="Size", xlab="Time")
}


# =============================================================================
# =============================================================================
# makeCovDf creates a dataframe of size variables for prediction
# TODO: make able to use 'covariates'
.makeCovDf <- function(size, explanatoryVariables) {
  sizeSorted <- sort(size)
  splitVars <- strsplit(explanatoryVariables, split = "\\+")
  expVar <- unlist(splitVars[grep("size", splitVars)])
  covDf <- data.frame(size = sizeSorted)
  for(i in 1:length(expVar)) {
    if(is.null(expVar[i])) next
    expVar[i] <- sub('[[:space:]]', '', expVar[i])
    if(expVar[i] == "size2") {
      covDf$size2 <- sizeSorted ^ 2
    }
    if(expVar[i] == "size3"){
      covDf$size3 <- sizeSorted ^ 3
    }
    if(expVar[i] == "expsize") {
      covDf$expsize <- exp(sizeSorted)
    }
    if(expVar[i] == "logsize") {
      covDf$logsize <- log(sizeSorted)
    }
    if(expVar[i] == "logsize2") {
      covDf$logsize2 <- log(sizeSorted ^ 2)
    }
  }
  return(covDf)
}

# =============================================================================
# =============================================================================
.makeListIPMs<- function(dataf, nBigMatrix=10, minSize=-2,maxSize=10, 
    integrateType="midpoint", correction="none",
    explSurv=surv~size+size2+covariates,
    explGrow=sizeNext~size+size2+covariates, 
    regType="constantVar",explFec=fec~size,Family="gaussian", 
    Transform="none",fecConstants=data.frame(NA)) {
  
  #convert to 1:n for indexing later
  dataf$covariates <- as.factor(dataf$covariates)
  levels(dataf$covariates) <- 1:length(unique(dataf$covariates))
  
  print(explSurv)
  sv1 <- makeSurvObj(dataf=dataf,
      Formula=explSurv)
  gr1 <- makeGrowthObj(dataf=dataf,
      Formula=explGrow,
      regType=regType)
  
  fv1 <- makeFecObj(dataf=dataf,Formula=explFec, Family=Family, Transform=Transform, 
      fecConstants=fecConstants) 
  
  covs <- unique(dataf$covariates)
  covs <- covs[!is.na(covs)]
  
  #print(covs)
  
  IPM.list <- list()
  for (k in 1:length(covs)) { 
    
    tpF <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
        maxSize = maxSize, chosenCov = data.frame(covariates=as.factor(k)),
        fecObj = fv1,integrateType=integrateType, correction=correction)
    tpS <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize,
        maxSize = maxSize, chosenCov = data.frame(covariates=as.factor(k)),
        growObj = gr1, survObj = sv1,
        integrateType=integrateType, correction=correction)
    IPM.list[[k]] <- tpF+tpS
  }
  return(IPM.list)
  
}

# =============================================================================
# =============================================================================
.convergeLambda<-function(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
    discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, tol=1e-4,binIncrease=5){
  
  lambda.new<-1000
  delta<-1000
  
  print(paste(c("delta: ","new lambda:", "New number of grid points:
                  "),collapse=""))
  
  
  while(delta>tol) {
    lambda.old <-lambda.new
    nBigMatrix <- nBigMatrix + binIncrease
    Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
        maxSize = maxSize, growObj = growObj, survObj = survObj, 
        discreteTrans = discreteTrans, integrateType = integrateType, 
        correction = correction)
    Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
        maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
        correction = correction, preCensus = preCensus,survObj=survObj,growObj=growObj)
    
    IPM <- Pmatrix + Fmatrix
    lambda.new <- Re(eigen(IPM)$value[1])
    
    delta<-abs(lambda.new-lambda.old)
    print(paste(c( delta, lambda.new,nBigMatrix),collapse=" "))
    
  }
  
  print(c("Final lambda from iteration:",lambda.new))
  print(c("Number of bins:",nBigMatrix))
  
  output <- list(nBigMatrix = nBigMatrix, IPM = IPM, lambda = lambda.new)
  
  return(output)
}


# =============================================================================
# =============================================================================
.convergeR0<-function(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
    discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, tol=1e-4,binIncrease=5){
  
  R0.new<-1000
  delta<-1000
  
  
  print(paste(c("delta: ","new R0:", "New number of grid points:"),collapse=""))
  
  while(delta>tol) {
    R0.old <-R0.new
    nBigMatrix <- nBigMatrix + binIncrease
    Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
        maxSize = maxSize, growObj = growObj, survObj = survObj, 
        discreteTrans = discreteTrans, integrateType = integrateType, 
        correction = correction)
    Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
        maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
        correction = correction, preCensus = preCensus,survObj=survObj,growObj=growObj)
    
    R0.new <- R0Calc(Pmatrix,Fmatrix)
    
    delta<-abs(R0.new-R0.old)
    print(paste(c( delta, R0.new,nBigMatrix),collapse=" "))
  }
  
  IPM <- Pmatrix + Fmatrix
  
  print(c("Final R0 from iteration:",R0.new))
  print(c("Number of bins:",nBigMatrix))
  
  output<-list(nBigMatrix = nBigMatrix, binIncrease=binIncrease,IPM=IPM,R0=R0.new)
  
  return(output)
}

# =============================================================================
# =============================================================================
.convergeLifeExpectancy <- function(growObj, survObj,nBigMatrix, minSize, maxSize, 
    discreteTrans = 1, integrateType = "midpoint", correction = "none", 
    tol=1e-1,binIncrease=5, chosenBin=1){
  
  LE.new<-1000
  delta<-1000
  
  print(paste(c("delta: ","new LE:", "New number of grid points:"),collapse=""))
  
  while(delta>tol) {
    LE.old <- LE.new
    nBigMatrix <- nBigMatrix + binIncrease
    Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
        maxSize = maxSize, growObj = growObj, survObj = survObj, 
        discreteTrans = discreteTrans, integrateType = integrateType, 
        correction = correction)
    
    LE.new <- meanLifeExpect(Pmatrix)    
    
    delta <- abs(LE.new[chosenBin]-LE.old[chosenBin])
    print(paste(c( delta, LE.new[chosenBin],nBigMatrix),collapse=" "))
  }
  
  print(c("Final life expectancy of chosen bin from iteration:",LE.new[1]))
  print(c("Number of bins:",nBigMatrix))
  
  output<-list(nBigMatrix = nBigMatrix,
      binIncrease=binIncrease,Pmatrix=Pmatrix,LE=LE.new)
  
  return(output)
}

# =============================================================================
# =============================================================================
# TO DO .fecPostCensus properly (i.e. include size changes between the census and reproduction event) the following steps have to be made:
# 1. use growSurv to determine what the distribution of size is at the reproduction event given initial size x
# 2. over this distribution of sizes x2 at the reproduction event, what are the expected number of offspring: per x2: .fecRaw(x=x2,...)
# 3. multiply the expected number of offspring per x2 with the probability that an offspring is of size y, using dnorm(y,predict(..., newdata=newd2, ...) where newd2 is calculated for all levels of x2
# all in all, Eelke wonders if the outer-solution is still useful in the .fecPostCensus case, since x2 needs to have a certain distribution, which is not passed down to .fecPostCensus... 

## REMOVE GROWTH FROM THIS - note that this means 
#### growth obj generally not needed down below.....
## A function that outer can use showing numbers from x to y via production, growth, survival and distribution offspring
.fecPostCensus <- function(x,y,cov=data.frame(covariate=1),
    fecObj, growObj,survObj,offspringObj=NULL) {
  newd <- data.frame(cbind(cov,size=x),
      stringsAsFactors = FALSE)
  
  newd$size2 <- x^2
  newd$size3 <- x^3
  
  if (is.null(offspringObj)){
    if (length(grep("expsize",fecObj@offspringRel$formula))>0 |
        length(grep("expsize",growObj@fit$formula))>0) { newd$expsize <- exp(x)}            	
    if (length(grep("logsize",fecObj@offspringRel$formula))>0 |
        length(grep("logsize",growObj@fit$formula))>0) { newd$logsize <- log(x)}            
    u <- .fecRaw(x=x,cov=cov,fecObj=fecObj)[[1]]*
        dnorm(y,predict(fecObj@offspringRel,newdata=newd, type="response"),fecObj@sdOffspringSize)*
        surv(size=x, cov=cov, survObj=survObj)
  } else {
    u <- .fecRaw(x=x,cov=cov,fecObj=fecObj)[[1]]*
        growth(y,y,newd,offspringObj)*surv(size=x, cov=cov, survObj=survObj)				
  }
  return(u)
}

# =============================================================================
# =============================================================================
## A function that outer can use giving pnorm for offspring reprod
.offspringCum <- function(x,y,cov=data.frame(covariate=1),fecObj,offspringObj=NULL) {
  newd <- data.frame(cbind(cov,size=x),
      stringsAsFactors = FALSE)
  
  if (is.null(offspringObj)){
    newd$size2 <- x^2
    newd$size3 <- x^3
    if (length(grep("expsize",fecObj@offspringRel$formula))>0) { newd$expsize <- exp(x)}            
    if (length(grep("logsize",fecObj@offspringRel$formula))>0) { newd$logsize <- log(x)}            
    u <- pnorm(y,predict(fecObj@offspringRel,newdata=newd, type="response"),fecObj@sdOffspringSize)
  } else { 
    u <- growthCum(y,y,newd, offspringObj)	
  }
  return(u)
}

# =============================================================================
# =============================================================================
#### growth obj generally not needed down below.....
## A function that outer can use showing numbers from x to y via production, growth, survival and distribution offspring
.fecPostCensusInteger <- function(x,y,cov=data.frame(covariate=1),fecObj, growObj,survObj, offspringObj=NULL) {
  newd <- data.frame(cbind(cov,size=x),
      stringsAsFactors = FALSE)
  
  u <- .fecRaw(x=x,cov=cov,fecObj=fecObj)[[1]]*
      surv(size=x, cov=cov, survObj=survObj)
  
  if (is.null(offspringObj)){
    newd$size2 <- x^2
    newd$size3 <- x^3
    if (length(grep("expsize",fecObj@offspringRel$formula))>0 |
        length(grep("expsize",growObj@fit$formula))>0) { newd$expsize <- exp(x)}            
    if (length(grep("logsize",fecObj@offspringRel$formula))>0 |
        length(grep("logsize",growObj@fit$formula))>0) { newd$logsize <- log(x)}            
    
    if (fecObj@distOffspring=="poisson")
      u <- u*dpois(y,predict(fecObj@offspringRel,newdata=newd, type="response"))
    if (fecObj@distOffspring=="negBin")
      u <- u*dnbinom(y,mu=predict(fecObj@offspringRel,newdata=newd, type="response"),
          size=fecObj@thetaOffspringSize)
  } else {
    u <- u*growth(y,y,newd, offspringObj)
    
  }
  
  
  return(u)
}
# =============================================================================
# =============================================================================

## A function that outer can use showing numbers from x to y via production and distribution offspring
.fecPreCensus <- function(x,y,cov=data.frame(covariate=1),fecObj,offspringObj=NULL) {
  newd <- data.frame(cbind(cov,size=x),
      stringsAsFactors = FALSE)	
  newd$size2 <- x^2
  newd$size3 <- x^3
  
  if (is.null(offspringObj)){
    if (length(grep("expsize",
            fecObj@offspringRel$formula))>0) { newd$expsize <- exp(x)}
    if (length(grep("logsize",
            fecObj@offspringRel$formula))>0) { newd$logsize <- log(x)}
    u <- .fecRaw(x=x,cov=cov,fecObj=fecObj)[[1]]*
        dnorm(y,predict(fecObj@offspringRel,newdata=newd, type="response"),
            fecObj@sdOffspringSize)
  } else {
    u <- .fecRaw(x=x,cov=cov,fecObj=fecObj)[[1]]*
        growth(y,y,newd,offspringObj)
    
  }
  #print(cbind(y,predict(fecObj@offspringRel)))
  
  
  return(u)
}

# =============================================================================
# =============================================================================
## A function that outer can use showing numbers from x to y via production and distribution offspring
.fecPreCensusInteger <- function(x,y,cov=data.frame(covariate=1),fecObj,offspringObj=NULL) {
  newd <- data.frame(cbind(cov,size=x),
      stringsAsFactors = FALSE)	
  newd$size2 <- x^2
  newd$size3 <- x^3
  
  u <- .fecRaw(x=x,cov=cov,fecObj=fecObj)[[1]]
  
  if (is.null(offspringObj)){
    if (length(grep("expsize",
            fecObj@offspringRel$formula))>0) { newd$expsize <- exp(x)}
    if (length(grep("logsize",
            fecObj@offspringRel$formula))>0) { newd$logsize <- log(x)}
    if (fecObj@distOffspring=="poisson")
      u <- u*dpois(y,predict(fecObj@offspringRel,newdata=newd, type="response"))
    if (fecObj@distOffspring=="negBin")
      u <- u*dnbinom(y,mu=predict(fecObj@offspringRel,newdata=newd, type="response"),
          size=fecObj@thetaOffspringSize)
  } else {
    u <- u*growth(y,y,newd, offspringObj)
    
  }
  
  #print(cbind(y,predict(fecObj@offspringRel)))
  
  return(u)
}

# =============================================================================
# =============================================================================
## Get raw numbers of offspring produced by every size class by multiplying up the constants,
## and doing all the "predict: values needed; and taking out only the babies that go to the continuous classes

.fecRaw <- function(x,cov=data.frame(covariate=1),fecObj) { 
  
  newd <- data.frame(cbind(cov,size=x),
      stringsAsFactors = FALSE)
  
  newd$size2 <- x^2
  newd$size3 <- x^3
  newd$expsize <- exp(x)
  
  #if (length(fecObj@offspringRel)>1) {
  #	if (length(grep("logsize", fecObj@offspringRel$formula))>0) { newd$logsize <- log(x)}
  #}
  
  fecObj@fecConstants[is.na(fecObj@fecConstants)] <- 1
  
  #fecundity rates
  fecValues <- matrix(c(rep(1,length(fecObj@fitFec)),unlist(fecObj@fecConstants)),
      ncol=length(x),nrow=length(fecObj@fitFec)+
          length(fecObj@fecConstants))
  #rownames(fecValues) <- c(fecObj@fecNames,names(fecObj@fecConstants))
  for (i in 1:length(fecObj@fitFec)) fecValues[i,] <- predict(fecObj@fitFec[[i]],newd,type="response")
  if (length(grep("log",fecObj@Transform))>0) for (i in grep("log",fecObj@Transform)) fecValues[i,]<-exp(fecValues[i,])
  if (length(grep("exp",fecObj@Transform))>0) for (i in grep("exp",fecObj@Transform)) fecValues[i,]<-log(fecValues[i,])
  if (length(grep("sqrt",fecObj@Transform))>0) for (i in grep("sqrt",fecObj@Transform)) fecValues[i,]<-(fecValues[i,])^2
  if (length(grep("-1",fecObj@Transform))>0) for (i in grep("-1",fecObj@Transform)) fecValues[i,]<-fecValues[i,]+1
  if (length(which(fecObj@vitalRatesPerOffspringType[,"continuous"]==1))>1) {
    prodFecValues <- apply(fecValues[which(fecObj@vitalRatesPerOffspringType[,"continuous"]==1),],2,prod)*unlist(fecObj@offspringSplitter["continuous"])
  } else {
    prodFecValues <- fecValues[which(fecObj@vitalRatesPerOffspringType[,"continuous"]==1),]*unlist(fecObj@offspringSplitter["continuous"])
  }
  return(list(prodFecValues,fecValues))
}

# =============================================================================
# =============================================================================

# function to predict growth object - NOTE THAT THIS will not work with factors / contrasts
.predictMuX <- function(grObj, newData, covPred = 1) {
  dataBase <- newData
  coefNames <- attr(grObj@fit$coefficients, "names")
  coefValues <- as.matrix(grObj@fit$coefficients)
  covNames <- coefNames[grep("covariate", coefNames)]
  covPos <- as.numeric(unlist(strsplit(covNames, "covariate")))
  covPos <- as.numeric(covPos[!is.na(covPos)])
  covDf <- as.data.frame(matrix(0, nrow = dim(newData)[1], ncol = length(covPos)))
  names(covDf) <- covNames
  # Turn "on" intended covariate and build newDataFrame
  if(covPred != 1) {
    covDf[, (covPred - 1)] <- 1
  }
  newData <- cbind(dataBase, covDf)
  newDataSubset <- as.matrix(cbind(1, newData[, (names(newData) %in% coefNames)]))
  predValues <- as.matrix(newDataSubset) %*% matrix(coefValues, length(coefValues), 1)
  return(as.numeric(predValues))
}

# =============================================================================
# =============================================================================
### sens params for discrete survival bit

.sensParamsDiscrete <-  function (growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
    chosenCov = data.frame(covariate=1),		
    discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, delta=1e-4) {
  
  
  #all the transitions - note that this is in the order of the columns, 
  nmes <- paste("",c(outer(colnames(discreteTrans@discreteTrans),colnames(discreteTrans@discreteTrans),paste,sep="-")),sep="")
  
  # all survival out of discrete stages
  nmes <- c(nmes,paste("survival from",colnames(discreteTrans@discreteSurv),sep=""))
  
  # all means out of discrete stages
  nmes <- c(nmes,paste("mean from ",colnames(discreteTrans@meanToCont),sep=""))
  
  # all sds out of discrete stages
  nmes <- c(nmes,paste("sd from ",colnames(discreteTrans@sdToCont),sep=""))
  
  #  not sure makes sense to do distribToDiscrete???	
  
  #coefficients linking continuous survival into discrete
  nmes <- c(nmes, paste("survival from continuous ", names(discreteTrans@moveToDiscrete$coefficients),sep=""))
  
  
  elam <- rep(NA,length(nmes))
  names(elam) <- nmes
  slam <- elam
  Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
      maxSize = maxSize, growObj = growObj, survObj = survObj, 
      chosenCov = chosenCov,
      discreteTrans = discreteTrans, integrateType = integrateType, 
      correction = correction)
  Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
      chosenCov = chosenCov,
      maxSize = maxSize, fecObj = fecObj, integrateType = integrateType, 
      correction = correction,preCensus = preCensus,survObj=survObj,growObj=growObj)
  IPM <- Pmatrix + Fmatrix
  lambda1 <- Re(eigen(IPM)$value[1])
  
  #all the transitions
  param.test <- 0
  for (j in 1:nrow(discreteTrans@discreteTrans)) {
    for (k in 1:nrow(discreteTrans@discreteTrans)) {
      
      #print(c(rownames(discreteTrans@discreteTrans)[k],colnames(discreteTrans@discreteTrans)[j]))
      
      param.test <- param.test+1
      
      if (discreteTrans@discreteTrans[k,j]==0) next()
      
      #pick element of matrix
      adj <- rep(discreteTrans@discreteTrans[k,j]*delta,nrow(discreteTrans@discreteTrans)-1)
      discreteTrans@discreteTrans[k,j] <- discreteTrans@discreteTrans[k,j] * (1 + delta)
      #alter the other values in the columns so as continue to sum to one
      if (sum(discreteTrans@discreteTrans[k,-j]>0)>0) adj <- adj/sum(discreteTrans@discreteTrans[k,-j]>0)
      adj[discreteTrans@discreteTrans[-k,j]==0] <- 0
      discreteTrans@discreteTrans[-k,j] <- discreteTrans@discreteTrans[-k,j]-adj
      
      #print(colSums(discreteTrans@discreteTrans))
      
      
      Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, 
          minSize = minSize, maxSize = maxSize, growObj = growObj, 
          chosenCov = chosenCov,
          survObj = survObj, discreteTrans = discreteTrans, 
          integrateType = integrateType, correction = correction)
      
      Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, 
          minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
          chosenCov = chosenCov,
          integrateType = integrateType, correction = correction,
          preCensus = preCensus,survObj=survObj,growObj=growObj)
      
      IPM <- Pmatrix + Fmatrix
      lambda2 <- Re(eigen(IPM)$value[1])
      discreteTrans@discreteTrans[k,j] <-  discreteTrans@discreteTrans[k,j]/(1 + delta)	
      discreteTrans@discreteTrans[-k,j] <- discreteTrans@discreteTrans[-k,j]+adj
      
      slam[param.test] <- (lambda2 - lambda1)/(discreteTrans@discreteTrans[k,j] * delta)
      elam[param.test] <- (lambda2 - lambda1)/(lambda1*delta)
    }}
  
  
  # all survival out of discrete stages
  count <- param.test
  for (param.test in 1:length(discreteTrans@discreteSurv)) { 
    
    discreteTrans@discreteSurv[param.test] <- discreteTrans@discreteSurv[param.test]* (1 + delta)
    
    Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, 
        minSize = minSize, maxSize = maxSize, growObj = growObj, 
        survObj = survObj, discreteTrans = discreteTrans, 
        chosenCov = chosenCov,
        integrateType = integrateType, correction = correction)
    
    Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, 
        minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
        chosenCov = chosenCov,
        integrateType = integrateType, correction = correction,
        preCensus = preCensus,survObj=survObj,growObj=growObj)
    
    IPM <- Pmatrix + Fmatrix
    lambda2 <- Re(eigen(IPM)$value[1])
    
    discreteTrans@discreteSurv[param.test] <- discreteTrans@discreteSurv[param.test]/ (1 + delta)
    slam[param.test+count] <- (lambda2 - lambda1)/(discreteTrans@discreteSurv[param.test] * delta)
    elam[param.test+count] <- (lambda2 - lambda1)/(lambda1*delta)
    
  }
  
  
  # all mean values coming out of discrete stages
  count <- param.test+count
  for (param.test in 1:length(discreteTrans@meanToCont)) { 
    
    discreteTrans@meanToCont[param.test] <- discreteTrans@meanToCont[param.test]* (1 + delta)
    
    Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, 
        minSize = minSize, maxSize = maxSize, growObj = growObj, 
        chosenCov = chosenCov,
        survObj = survObj, discreteTrans = discreteTrans, 
        integrateType = integrateType, correction = correction)
    
    Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, 
        minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
        chosenCov = chosenCov,
        integrateType = integrateType, correction = correction,
        preCensus = preCensus,survObj=survObj,growObj=growObj)
    
    IPM <- Pmatrix + Fmatrix
    lambda2 <- Re(eigen(IPM)$value[1])
    
    discreteTrans@meanToCont[param.test] <- discreteTrans@meanToCont[param.test]/ (1 + delta)
    slam[param.test+count] <- (lambda2 - lambda1)/(discreteTrans@meanToCont[param.test] * delta)
    elam[param.test+count] <- (lambda2 - lambda1)/(lambda1*delta)
    
  }
  
  
  # all sd values coming out of discrete stages
  count <- param.test+count
  for (param.test in 1:length(discreteTrans@sdToCont)) { 
    
    discreteTrans@sdToCont[param.test] <- discreteTrans@sdToCont[param.test]* (1 + delta)
    
    Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, 
        minSize = minSize, maxSize = maxSize, growObj = growObj, 
        chosenCov = chosenCov,
        survObj = survObj, discreteTrans = discreteTrans, 
        integrateType = integrateType, correction = correction)
    
    Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, 
        minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
        chosenCov = chosenCov,
        integrateType = integrateType, correction = correction,
        preCensus = preCensus,survObj=survObj,growObj=growObj)
    
    IPM <- Pmatrix + Fmatrix
    lambda2 <- Re(eigen(IPM)$value[1])
    
    discreteTrans@sdToCont[param.test] <- discreteTrans@sdToCont[param.test]/ (1 + delta)
    slam[param.test+count] <- (lambda2 - lambda1)/(discreteTrans@sdToCont[param.test][param.test] * delta)
    elam[param.test+count] <- (lambda2 - lambda1)/(lambda1*delta)
    
  }
  
  # parameters linking size to survival
  count <- param.test+count
  for (param.test in 1:length(discreteTrans@moveToDiscrete$coefficients)) { 
    
    discreteTrans@moveToDiscrete$coefficients[param.test] <- discreteTrans@moveToDiscrete$coefficients[param.test]* (1 + delta)
    
    Pmatrix <- makeIPMPmatrix(nBigMatrix = nBigMatrix, 
        minSize = minSize, maxSize = maxSize, growObj = growObj, 
        chosenCov = chosenCov,
        survObj = survObj, discreteTrans = discreteTrans, 
        integrateType = integrateType, correction = correction)
    
    Fmatrix <- makeIPMFmatrix(nBigMatrix = nBigMatrix, 
        minSize = minSize, maxSize = maxSize, fecObj = fecObj, 
        chosenCov = chosenCov,
        integrateType = integrateType, correction = correction,
        preCensus = preCensus,survObj=survObj,growObj=growObj)
    
    IPM <- Pmatrix + Fmatrix
    lambda2 <- Re(eigen(IPM)$value[1])
    
    discreteTrans@moveToDiscrete$coefficients[param.test] <- discreteTrans@moveToDiscrete$coefficients[param.test]/ (1 + delta)
    slam[param.test+count] <- (lambda2 - lambda1)/(discreteTrans@moveToDiscrete$coefficient[param.test] * delta)
    elam[param.test+count] <- (lambda2 - lambda1)/(lambda1*delta)
    
  }
  
  print("Did not calculate sensitivities and elasticities for:")
  print(c(names(slam[is.na(slam)])))
  print("Values of zero")
  
  slam <- slam[!is.na(slam)]
  elam <- elam[!is.na(elam)]
  
  return(list(slam = slam, elam = elam))
  
}

# =============================================================================
# =============================================================================
#Generic for life expectancy in a modeled stoch env  -NEEDS DOUBLE-CHECKING
#parameters - a compound IPM of dim nEnvClass*nBigMatrix
#           - an environmental matrix
# returns - the life expectancy for each of the sizes in the IPM (columns)
#           for each of the starting env states
#
#
.stochLifeExpect <- function(IPMmatrix,envMatrix){
  
  matrix.dim <- length(IPMmatrix[1,])
  nstages <- IPMmatrix@nBigMatrix
  nstates <- IPMmatrix@nEnvClass
  
  pis <- Re(eigen(envMatrix)$vector[,1])
  pis <- pis/(sum(pis))
  
  #print(pis)
  
  #ckron <- kronecker(envMatrix,diag(nstages))  #doesnt work??
  m <- IPMmatrix  #%*%ckron  #eq 29 in carols paper
  
  Itilda <- diag(matrix.dim)
  
  #not used in this one
  eatildas <- array(dim=c(nstates*nstages,nstages,nstates))
  eatildas[,,] <-0
  for (i in 1:nstates){
    eatildas[,,i] <- Itilda[,((i-1)*nstages+1):(nstages*i)]
  }
  
  qatildas <- array(dim=c(nstates*nstages,nstages,nstates));
  qatildas[,,]<-0
  for (i in 1:nstates) {
    indext <-((i-1)*nstages+1):(nstages*i) 
    qatildas[indext,,i] <- IPMmatrix[indext,indext]/envMatrix[i,i]
    #print( IPMmatrix[cbind(indext,indext)]/envMatrix[i,i] )
  }                            #array of qatildas, eqn 26
  #need to remove env effect since mega-matrix pre-built 
  
  #print(qatildas)
  
  etilda <- array(dim=c(nstates*nstages,nstages));
  etilda[,]<-0
  for (i in 1:nstates) { 
    etilda[((i-1)*nstages+1):(nstages*i),] <- diag(nstages);
  }                            #etilda, eqn 27
  
  I <- diag(nstages);                 #identity matrix of dimension K x K
  
  Ns.markov <- array(dim=c(nstages,nstages,nstates)); #set up for array of Ns
  Ns.markov[,,]<-0
  for (i in 1:nstates){ 
    Ns.markov[,,i] <- I + t(etilda)%*%(solve(Itilda-m))%*%qatildas[,,i];
    #eqn 28, conditional on initial state 1
  }
  
  #print(Ns.markov)
  
  #average over all initial states %%%%%
  Nbar.markov <- rep(0,nstages);
  for (i in 1:nstates) { 
    Nbar.markov <-  Nbar.markov + pis[i] * Ns.markov[,,i]; 
  }                         #eqn 29, weight each fundamntal matrix by expctd frequency 
  
  
  
  #lifeexp, column sums
  lifeexp.markov<-matrix(0,nstates,nstages); #set up array
  for (i in 1:nstates) { 
    lifeexp.markov[i,] <- apply(Ns.markov[,,i], 2, sum);
  }
  
  return(lifeexp.markov)
}

# =============================================================================
# =============================================================================
# =============================================================================
# =============================================================================
# Function to augment mortality data for large trees
# Parameters dataf - existing data frame
#            size.thresh - the size above which tree death is being augmented
#            prop.dead - the proportion of these expected to be dead

.deathDataAugment <- function (dataf, size.thresh, prop.dead) { 
  
  n.now <- sum(dataf$size > size.thresh)
  n.new.dead <- ceiling(prop.dead * n.now / (1 - prop.dead))
  new.size <- rnorm(n.new.dead,size.thresh, sd(dataf$size[dataf$size > size.thresh]))
  
  datanew <- data.frame(size = new.size, sizeNext = rep(NA, n.new.dead), surv = rep(0, n.new.dead), 
      covariate = rep(0, n.new.dead), covariateNext = rep(0, n.new.dead),
      fec = rep(NA, n.new.dead)) 
  
  dataf.new <- rbind(dataf,datanew)
  
  return(dataf.new)
  
}

# =============================================================================
# =============================================================================
# replace the growth object fit with a new, desired variance for predict
.alteredFit <- function(dummyFit = dummyFit, 
    newCoef = dummyFit$coefficients, 
    desiredSd = 1) {
  dummyFit$coefficients[] <- newCoef
  dummyFit$residuals <- rnorm(length(dummyFit$residuals), mean = 0, sd = desiredSd)	
  # need to use qr here to assign dummyFit so that there is no warning when decomposed for n
  # Error in rnorm(residDf, 0, sd = desiredSd) : object 'residDf' not found
  return(dummyFit)	
}

# =============================================================================
# =============================================================================
# Function to take a list of growth and survival objects and make a list of Pmatrices
#
# Parameters - growObjList - a list of growth objects
#            - survObjList - a list of survival objects
#            - nBigMatrix - the number of bins
#            - minSize - the minimum size
#            - maxSize - the maximum size
#            - cov - is a discrete covariate considered
#            - envMat - enviromental matrix for transition between
# 
# Returns    - a list of Pmatrices
.makeListPmatrix  <- function(growObjList,survObjList,
    nBigMatrix,minSize,maxSize, cov=FALSE, envMat=NULL,
    integrateType="midpoint",correction="none", discreteTransList=1) {
  
  if (length(growObjList)>length(survObjList)) { 
    survObjList <- sample(survObjList,size=length(growObjList),replace=TRUE)
  } else { 
    if (length(growObjList)<length(survObjList))  
      growObjList <- sample(growObjList,size=length(survObjList),replace=TRUE)
  }
  
  nsamp <- length(growObjList)
  
  if(length(discreteTransList)<nsamp ){
    # if(warn) warning('Length of discreteTrans list is less than the length of another vital rate object list, so some members of the discreteTrans list have been repeated.')
    discreteTransList <- sample(discreteTransList,size=nsamp,replace=T)
  }
  
  PmatrixList <- list()
  for ( k in 1:length(growObjList)) { 
    if (!cov) {
      PmatrixList[[k]] <- makeIPMPmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
          maxSize = maxSize, growObj = growObjList[[k]],discreteTrans=discreteTransList[[k]],
          survObj = survObjList[[k]],integrateType=integrateType, correction=correction) 
    } else {
      PmatrixList[[k]] <- makeCompoundPmatrix(nEnvClass = length(envMat[1,]),
          nBigMatrix = nBigMatrix, minSize = minSize, 
          maxSize = maxSize, envMatrix=envMat,discreteTrans=discreteTransList[[k]],
          growObj = growObjList[[k]],
          survObj = survObjList[[k]],integrateType=integrateType, correction=correction)    
    }
  }
  
  return(PmatrixList)
}

# =============================================================================
# =============================================================================
# Function to take a list of growth and survival objects and make a list of Fmatrices

.makeListFmatrix  <- function(fecObjList,nBigMatrix,minSize,maxSize, cov=FALSE, 
    envMat=NULL,integrateType="midpoint",correction="none") {
  
  nsamp <- max(length(fecObjList))
  if (length(fecObjList)<nsamp)  
    fecObjList <- sample(fecObjList,size=nsamp,replace=TRUE)
  
  FmatrixList <- list()
  for ( k in 1:nsamp) {
    if (!cov) { 
      FmatrixList[[k]] <- makeIPMFmatrix(nBigMatrix = nBigMatrix, minSize = minSize, 
          maxSize = maxSize, 
          fecObj=fecObjList[[k]],integrateType=integrateType, correction=correction)
    } else {
      FmatrixList[[k]] <- makeCompoundFmatrix(nEnvClass = length(envMat[1,]),
          nBigMatrix = nBigMatrix, minSize = minSize, 
          maxSize = maxSize, envMatrix=envMat,
          fecObj=fecObjList[[k]],integrateType=integrateType, correction=correction)
    }
    
    FmatrixList[[k]] <-  FmatrixList[[k]]
  }
  return(FmatrixList)
}


# =============================================================================
# =============================================================================
### new functions createGrowhtObj and createSurvObj which will make up their own data

.createGrowthObj <- function(Formula=sizeNext~size, coeff=c(1,1), sd=1){ 
  
  var.names <- all.vars(Formula)
  
  if (length(coeff)!=(length(var.names))) #length var.names, because intercept not named 
    stop("not enough coefficients supplied for the chosen Formula")
  
  dataf<- as.data.frame(matrix(rnorm(3*length(var.names)),3,length(var.names)))
  colnames(dataf) <- var.names
  
  fit <- lm(Formula, data=dataf)
  
  if (length(grep("sizeNext", as.character(Formula))) > 0) { 
    
    gr1 <- new("growthObj")
    gr1@fit <- fit
    gr1@fit$coefficients <- coeff
    gr1@sd <- sd
  }  
  
  if (length(grep("incr", as.character(Formula))) > 0) { 
    
    gr1 <- new("growthObjIncr")
    gr1@fit <- fit
    gr1@fit$coefficients <- coeff
    gr1@sd <- sd
  }  
  
  return(gr1)
  
}

# =============================================================================
# =============================================================================
.createSurvObj <- function(Formula=surv~size, coeff=c(1,1)){ 
  var.names <- all.vars(Formula)
  
  #not that although var.names will have one extra (cos of response variable
  # this will correspond to the intercept
  if (length(coeff)!=(length(var.names))) 
    stop("not enough coefficients supplied for the chosen Formula")
  
  dataf<- as.data.frame(matrix(rnorm(3*length(var.names)),3,length(var.names)))
  colnames(dataf) <- var.names
  dataf$surv <- sample(c(0,1),nrow(dataf), replace=TRUE)
  
  fit <- glm(Formula, data=dataf, family=binomial)
  
  sv1 <- new("survObj")
  sv1@fit <- fit
  
  return(sv1)
  
}

# =============================================================================
# =============================================================================
.createFecObj <- function(Formula=list(fec1~size,fec2~size+size2), 
    coeff=list(c(1,1),c(1,1,1)),
    Family = c("gaussian","binomial"),
    Transform = c("log","none"),
    meanOffspringSize = NA, sdOffspringSize = NA, 
    offspringSplitter = data.frame(continuous = 1), 
    vitalRatesPerOffspringType = data.frame(NA), 
    fecByDiscrete = data.frame(NA), 
    offspringSizeExplanatoryVariables = "1",
    fecConstants = data.frame(NA), 
    doOffspring=TRUE, 
    reproductionType="sexual"){ 
  var.names <- c()
  fecNames <- rep(NA,length(Formula))
  
  for (j in 1:length(Formula)) { 
    
    
    fecNames[j] <- all.vars(Formula[[j]])[1]
    var.names.here <- all.vars(Formula[[j]])
    
    if (length(coeff[[j]])!=(length(var.names.here))) 
      stop(paste("not enough coefficients supplied for the ",j, "th Formula", sep=""))
    
    var.names <- c(var.names,var.names.here)
  }
  
  var.names <- unique(var.names)
  
  #build a data-frame with all the right variables
  dataf<- as.data.frame(matrix(rnorm(3*length(var.names)),3,length(var.names)))
  colnames(dataf) <- var.names
  dataf$surv <- sample(c(0,1),nrow(dataf), replace=TRUE)
  
  if (sum(Transform=="log")>0) dataf[,fecNames[which(Transform=="log")]] <- pmax(dataf[,fecNames[which(Transform=="log")]],1)
  if (sum(Family=="binomial")>0) dataf[,fecNames[which(Family=="binomial")]] <- rbinom(nrow(dataf),1,0.5)
  
  fv1 <- makeFecObj(dataf=dataf, 
      fecConstants = fecConstants, 
      Formula = Formula, 
      Family = Family, 
      Transform = Transform, 
      meanOffspringSize = meanOffspringSize, 
      sdOffspringSize = sdOffspringSize, offspringSplitter = offspringSplitter, 
      vitalRatesPerOffspringType = vitalRatesPerOffspringType, 
      fecByDiscrete = fecByDiscrete, 
      offspringSizeExplanatoryVariables = offspringSizeExplanatoryVariables, 
      coeff=NULL, doOffspring=doOffspring, 
      reproductionType=reproductionType)
  
  #now over-write with the desired coefficients!
  for (j in 1:length(Formula)) { 
    fv1@fitFec[[j]]$coefficients <- coeff[[j]]
  }
  
  return(fv1)
  
}

# =============================================================================
# =============================================================================






















