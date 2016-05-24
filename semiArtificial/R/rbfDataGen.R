
rbfDataGen <- function(formula, data, eps=1e-4, minSupport=1, nominal=c("encodeBinary", "asInteger")) {
  
  if (!inherits(formula,"formula")) 
    stop("First argument must be a formula.");
  
  #minSupport = max(minSupport, 2)
  nominal <- match.arg(nominal)
  
  #prepare data according to the supplied formula, dependent variable on place 1
  dat <- model.frame(formula, data=data, na.action=na.pass);
  
  datNames <- names(dat)
  originalNames <- names(data)[sort(match(datNames, names(data)))]
  noInst <- nrow(dat)
  noAttr = ncol(dat) - 1
  
  attrClasses<-list()
  attrLevels <-list()
  attrOrdered<-logical()
  noCol <-vector(mode="integer",length=ncol(dat)) # number of columns when treating nominal attributes "encodeBinary"
  isDiscrete <- vector(mode="logical",length=ncol(dat))
  for (i in 1:ncol(dat)) {
    attrClasses[[i]] <- class(dat[[i]])
    attrOrdered[i] <- is.ordered(dat[[i]])
    if (is.factor(dat[[i]])){
      isDiscrete[i] <- TRUE
      attrLevels[[i]] <- levels(dat[[i]])
      noCol[i] <- length(attrLevels[[i]])
      if (noCol[i] == 2)
        noCol[i] <- 1
    } 
    else {
      attrLevels[[i]] <- NULL    
      noCol[i] <- 1
    }
  }
  noCol[1] <- 1 # class remains unchanged, even if it is a factor with several values
  noClasses <- length(attrLevels[[1]])
  classProb <- table(dat[[1]])/nrow(dat)
     
  if (nominal=="encodeBinary") {   
    isDiscrete <- vector(mode="logical", length= sum(noCol))
    datBin <- as.data.frame(matrix(0,nrow=nrow(dat),ncol=sum(noCol)))
    datBin[[1]] <- dat[[1]] # class is not transformed
    isDiscrete[1] <- TRUE
    idx <- 2    
    for (i in 2:ncol(dat)) {
      if (is(dat[[i]], "factor")) {
        isDiscrete[idx:(idx+noCol[i]-1)] <- TRUE
        if (noCol[i] > 1) {
           datBin[,idx:(idx+noCol[i]-1)] <- decodeClassLabels(dat[,i])
        }
        else datBin[[idx]] <- as.integer(dat[,i])
      }
      else {
        datBin[[idx]] <- dat[[i]]
      }
      idx <- idx + noCol[i]  
    }
    dat <- datBin
  }
  else if (nominal=="asInteger") {
    for (i in 2:ncol(dat)) {
      if (is(dat[[i]], "factor")) 
        dat[[i]] <- as.integer(dat[[i]])
      }
  }

  noAttrGen <- ncol(dat)-1
  # shuffle the data
  dat <- dat[sample(1:nrow(dat),nrow(dat)),]
  # prepare data for RSNNS model
  dataLearn <- prepareRSNNS(dat, targetCol=1, ratio=0) 
  # learn RBF model 
  model <- rbfDDA(dataLearn$inputsTrain, dataLearn$targetsTrain, linOut=TRUE)
  
  # extract network structure and parameters
  hidden <- model$snnsObject$getAllHiddenUnits()
  bias <- model$snnsObject$getUnitDefinitions()[model$snnsObject$getAllHiddenUnits(),"unitBias"]
  centers <- model$snnsObject$getWeightMatrix(model$snnsObject$getAllInputUnits(),model$snnsObject$getAllHiddenUnits())
  weights <- apply(model$snnsObject$getWeightMatrix(model$snnsObject$getAllHiddenUnits(),model$snnsObject$getAllOutputUnits()),1,sum)
  if (max(weights) < minSupport) {
    minSupport <- 1
    warning("rbfDataGen: all kernel weights are bellow minSupport, setting minSupport to 1 ")
  }
  probs <- weights
  probs[weights < minSupport] <- 0
  probs <- probs / sum(probs)
  unitClass <- apply(model$snnsObject$getWeightMatrix(model$snnsObject$getAllHiddenUnits(),model$snnsObject$getAllOutputUnits()),1,which.max)
  
  noG <- length(hidden)
  #assign each training instance to one cluster and compute sigma based on that
  assignedInst <- list()
  for (g in 1:noG) {  
    assignedInst[[g]] <- vector(mode="numeric", length=0)
  }
  activation <- vector(mode="numeric", length=noG)
  for (i in 1:noInst) {
    for (j in 1:noG) {
      activation[j] <- exp(-sum((dataLearn$inputsTrain[i,]-centers[,j])^2)*bias[j])
    }
    # assign to Gausian with higher activation level
    maxAct <- which.max(activation)
    assignedInst[[maxAct]]<-c(assignedInst[[maxAct]],i)
  }
  # compute spread of assigned instances around center
  spread <- matrix(0, nrow=noG,ncol=noAttrGen)
  gNoActivated = c()
  for (g in 1:noG) {  
    gNoActivated[g] <- length(assignedInst[[g]])
    for (a in 1: noAttrGen){
      spread[g, a] <- cspread(centers[a,g], dataLearn$inputsTrain[assignedInst[[g]],a]) 
      if (spread[g,a] < eps || is.nan( spread[g, a]) || is.infinite(spread[g, a]) )
        spread[g,a] <- eps
    }
  }
  gen <- list(noGaussians=noG, centers=centers, probs=probs, unitClass=unitClass, bias=bias,
              spread=spread, gNoActivated=gNoActivated, noClasses=noClasses, classProb = classProb,
              noAttr=noAttr, datNames=datNames,  originalNames=originalNames, 
              attrClasses=attrClasses, attrLevels=attrLevels, attrOrdered=attrOrdered,
              normParameters=getNormParameters(dataLearn$inputsTrain),
              noCol=noCol,isDiscrete=isDiscrete, noAttrGen=noAttrGen, nominal=nominal,
              eps=eps)
  class(gen) <- "RBFgenerator"
  gen
}

newdata <- function(object, ...) UseMethod("newdata", object)

newdata.RBFgenerator <- function(object, size, var=c("estimated","Silverman"), classProb=NULL, defaultSpread=0.05, ...) {
  if (! is(object, "RBFgenerator"))
    stop("newdata requires object of type RBFgenerator")
  
  var <- match.arg(var)
  n <- size
  
  if (is.null(classProb)) { # use class probability from the generator
    classProb <- object$classProb
  }
  
  # compute number of cases to generate from each Gaussian
  if (length(classProb) != object$noClasses)
      warning("Class probability distribution does not match the number of classes in the generator which is ", object$noClasses)
  if (abs(sum(classProb)-1) > object$eps)
        warning("Class probability distribution does not sum to 1 for classProb ")
  unitProbs <- vector(mode="numeric", length = object$noGaussians)
  for (cl in 1:object$noClasses) {
       sumProb <- sum(object$probs[object$unitClass == cl])
       if (sumProb > object$eps)
          unitProbs[object$unitClass == cl] <- classProb[cl] * object$probs[object$unitClass == cl] / sumProb             
  }
  # some classes might not be represented with the kernels, so we renormalize the distribution
  unitProbs <- unitProbs / sum(unitProbs)
  
  nh <- intFromProb(unitProbs, n) # number of instances from each Gaussian unit
  
  Rprof(memory.profiling=TRUE)
  genValues <- matrix(0,nrow=0,ncol=object$noAttrGen)
  genClass <- c()
  sigma <- diag(nrow=object$noAttrGen)
  
  spread <- object$spread
  if (!is.null(defaultSpread))
      spread[spread <= object$eps] <- defaultSpread
  
  for (i in 1:object$noGaussians) {
    if (nh[i] > 0) {
      # prepare  diagonal matrix Sigma
      if (var=="Silverman") { # Silverman's rule of a thumb 
        #  spread[g,a] contains variance of values in dimension a for kernel g
        #for (a in 1: object$noAttrGen)
        #  spread[i, a] = (4/(object$noAttrGen+2))^(1/(object$noAttrGen+4)) * object$gNoActivated[i]^(-1/(object$noAttrGen+4)) * sqrt(spread[i,a])
        spread[i,] = (4/(object$noAttrGen+2))^(1/(object$noAttrGen+4)) * object$gNoActivated[i]^(-1/(object$noAttrGen+4)) * sqrt(spread[i,])
      }
      # for "estimated" spread of training instances serves as variance, for "Silverman" his rule multivariate rule of a thumb is used
      diag(sigma) <- spread[i,]
      multiplier <- 2
      ggenerated <- 0
      while ( ggenerated < nh[i] ) {
        gvi <- mvrnorm(n=max(10,multiplier*nh[i]), mu=object$centers[,i], Sigma=sigma) # run generator
        gvi[,object$isDiscrete[-1]] <- round(gvi[,object$isDiscrete[-1]]) # round all discrete attributes to integers
        #activation <- c()
        #for (case in 1: nrow(gvi)) 
        #  activation[case] <- exp(-bias[i]*(gvi[case,]-centers[,i])^2)
        #gvi <- gvi[activation>=0.4,,drop=F] # remove instances with activation below 0.4
        rowmin <- apply(gvi, 1, min) # min row values
        rowmax <- apply(gvi, 1, max) # max row values
        # remove instances out of normalized rang; for high dimensional data sets this might significantly reduce
        # the number of acceptable instances, therefore we repeat this in a loop
        if (multiplier <= 128)
          gvi <- gvi[rowmin >=0 & rowmax <= 1, , drop=F]
        else {
          gvi[gvi > 1] <- 1
          gvi[gvi < 0] <- 0       
        }
        if (nrow(gvi) > nh[i]-ggenerated)
          gvi <- gvi[1:(nh[i]-ggenerated), ,drop=F]
        ggenerated <- ggenerated + nrow(gvi)
        genValues <- rbind(genValues, gvi)
        genClass <- c(genClass,rep(object$unitClass[i],nrow(gvi)))
        multiplier <- 2 * multiplier
        if (multiplier > 128) {
          warning("newdata cannot generate sufficient data from Gaussian number ", i)
        }
      }
    }
  }
  # denormalize
  genValues <- denormalizeData(genValues, object$normParameters)
  if (object$nominal=="encodeBinary") {
    genAttr <- matrix(0,nrow=nrow(genValues),ncol=object$noAttr)
    idx <- 1
    for (i in 1:ncol(genAttr)) {
      if ( object$noCol[i+1] > 1 )  
         genAttr[,i] <- binary2nominal(genValues[,idx:(idx+object$noCol[i+1]-1)]) 
      else genAttr[,i] <- genValues[,idx]
      idx <- idx + object$noCol[i+1]                 
    }
    genValues <- genAttr
  }
  
  genData<-as.data.frame(matrix(cbind(genClass,genValues),nrow=length(genClass),ncol=ncol(genValues)+1))
 
  for (i in 1:ncol(genData)){
     if ( "factor" %in% object$attrClasses[[i]])  {
        genData[[i]] <- factor(genData[[i]],levels=1:length(object$attrLevels[[i]]), labels=object$attrLevels[[i]])
        if (object$attrOrdered[i])
          genData[[i]] <- as.ordered(genData[[i]])
    }
  }
  names(genData) <- object$datNames
  genData <- genData[sample(1:nrow(genData),nrow(genData)),]  # shuffle
  Rprof(NULL)
  genData[,object$originalNames]
}  
  

#prepare format for RSNNS
prepareRSNNS <- function(data, targetCol=ncol(data),ratio=0.3){
  data <- medianImpute(data) # impute missing values with column median or mode 
  dataValues <- data[,-targetCol,drop=F]
  dataTargets <- decodeClassLabels(data[,targetCol])  
  dataLearn <- splitForTrainingAndTest(dataValues, dataTargets, ratio=ratio) # ratio% of test data
  dataLearn <- normTrainingAndTestSet(dataLearn,type="0_1") # normalize to  interval [0, 1]
  dataLearn
}

cspread <- function(center, values) {
  sum((values-center)^2)/(length(values)-1)
}

nominal2binary <- function(attr) {
  levels <- levels(attr)
  noVals <- length(levels)
  if (length(levels) <= 2)
    return(attr)
  else {
    bindata <- as.data.frame(matrix(FALSE, nrow <- length(attr), ncol= noVals))
    for (i in 1:ncol(bindata)) {
      bindata[,i] <- as.integer(attr)==i    
      bindata[[i]] <- factor(as.integer(bindata[[i]]),levels=0:1)
    }
    return(bindata)
  }
}

binary2nominal <- function(bindata) {
  attr <- vector(mode="integer",length=nrow(bindata))
  attr <- apply(bindata,1,which.is.max)

}

