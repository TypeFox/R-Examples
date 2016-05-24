# TODO: Add comment
# 
###############################################################################


## FUNCTIONS TO BE DEPRECATED #################################################
# =============================================================================
# =============================================================================
createIPMFmatrix <- function(fecObj,
    nEnvClass = 1,
    nBigMatrix = 50,
    minSize = -1,
    maxSize = 50,
    chosenCov = data.frame(covariate=1),
    integrateType="midpoint",
    correction="none",
    preCensus=TRUE,
    survObj=NULL,
    growObj=NULL) {
  
  print("warning: deprecated! use makeIPMFmatrix instead")
  
  rc <- makeIPMFmatrix(fecObj=fecObj,
      nEnvClass = nEnvClass,
      nBigMatrix = nBigMatrix,
      minSize = minSize,
      maxSize = maxSize,
      chosenCov = chosenCov,
      integrateType=integrateType,
      correction=correction,
      preCensus=preCensus,
      survObj=survObj,
      growObj=growObj)
  
  return(rc)
  
}

# =============================================================================
# =============================================================================
createIPMCmatrix <- function(clonalObj,
    nEnvClass = 1,
    nBigMatrix = 50,
    minSize = -1,
    maxSize = 50,
    chosenCov = data.frame(covariate=1),
    integrateType="midpoint",
    correction="none",
    preCensus=TRUE,
    survObj=NULL,
    growObj=NULL) {
  
  print("warning: deprecated! use makeIPMCmatrix instead")
  
  rc <- makeIPMFmatrix(fecObj=clonalObj,
      nEnvClass = nEnvClass,
      nBigMatrix = nBigMatrix,
      minSize = minSize,
      maxSize = maxSize,
      chosenCov = chosenCov,
      integrateType=integrateType,
      correction=correction,
      preCensus=preCensus,
      survObj=survObj,
      growObj=growObj)
  
  return(rc)
  
}

# =============================================================================
# =============================================================================
createIPMPmatrix <- function (nEnvClass = 1, 
    nBigMatrix = 50, minSize = -1, maxSize = 50, 
    chosenCov = data.frame(covariate = 1), growObj, survObj, 
    discreteTrans = 1, integrateType = "midpoint", correction = "none") {
  
  print("warning: deprecated! use makeIPMPmatrix instead")
  
  rc <- makeIPMPmatrix(nEnvClass = nEnvClass,
      nBigMatrix = nBigMatrix,
      minSize = minSize,
      maxSize = maxSize,
      chosenCov = chosenCov,
      growObj=growObj, survObj=survObj,
      discreteTrans =discreteTrans,
      integrateType=integrateType,
      correction=correction)
  
  return(rc)
}


# =============================================================================
# =============================================================================
createIntegerPmatrix <- function (nEnvClass = 1, 
    meshpoints=1:20,
    chosenCov = data.frame(covariate = 1), 
    growObj, survObj, 
    discreteTrans = 1){
  
  print("warning: deprecated! use makeIntegerPmatrix instead")
  
  rc <- makeIntegerPmatrix(nEnvClass = nEnvClass,meshpoints=meshpoints,
      chosenCov = chosenCov,growObj=growObj, survObj=survObj,discreteTrans =discreteTrans)
  
  return(rc)
  
}

# =============================================================================
# =============================================================================
createIntegerFmatrix <- function(fecObj,
    nEnvClass = 1,
    meshpoints=1:20,
    chosenCov = data.frame(covariate=1),
    preCensus=TRUE,
    survObj=NULL,
    growObj=NULL){
  
  print("warning: deprecated! use makeIntegerFmatrix instead")
  
  rc <- makeIntegerFmatrix(fecObj=fecObj,
      nEnvClass = nEnvClass,meshpoints=meshpoints,
      chosenCov = chosenCov,preCensus=preCensus, 
      growObj=growObj, survObj=survObj)
  
  
  return(rc)
  
}

# =============================================================================
# =============================================================================
#Function creates a combo of F.IPMs for a chosen
#size range, with range env categories and transitions 
#between them
#and growth and survival objects
#
#Parameters -
#   nEnvClass - the number of env classes, defaults to 2, should match dim envMatrix
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   envMatrix - a matrix describing transiions between env
#   fecObj - a fecundity object
#   integrateType - NOT YET IMPLEMENTED
#
#Returns -
#  an IPM object


makeCompoundCmatrix <- function(nEnvClass = 2,
    nBigMatrix = 50,
    minSize = -1,
    maxSize = 50,
    envMatrix,
    clonalObj,
    integrateType="midpoint",
    correction="none",
    preCensus=TRUE,
    survObj=NULL,
    growObj=NULL, offspringObj= NULL) {
  
  rc <- makeCompoundFmatrix(nEnvClass = nEnvClass,
      nBigMatrix = nBigMatrix,
      minSize = minSize,
      maxSize = maxSize,
      envMatrix = envMatrix,
      fecObj = clonalObj,
      integrateType=integrateType,
      correction=correction,
      preCensus=preCensus,
      survObj=survObj,
      growObj=growObj, offspringObj= offspringObj)
  
  return(rc) 
}

# =============================================================================
# =============================================================================
### CREATE DISCRETE F MATRIX
makeIntegerFmatrix <- function(fecObj,
    nEnvClass = 1,
    meshpoints=1:20,
    chosenCov = data.frame(covariate=1),
    preCensus=TRUE,
    survObj=NULL,
    growObj=NULL, offspringObj=NULL) {
  
  # boundary points b and mesh points y
  y <- meshpoints;
  nBigMatrix <- length(y)
  
  fecObj@fecConstants[is.na(fecObj@fecConstants)] <- 1
  
  tmp <- matrix(0,ncol=length(y),nrow=length(y))
  
  # 1. pre-census
  if (preCensus) { 
    ##NOTE that the condition is necessary cos you might ONLY have discrete offspring
    if (fecObj@offspringSplitter$continuous>0) { 
      tmp <- t(outer(X=y,Y=y,.fecPreCensusInteger,cov=chosenCov,fecObj=fecObj, offspringObj= offspringObj))
    }
    ##NOTE that the condition is necessary because you might ONLY have discrete offspring
    
    
    # 2. post-census
  } else {
    # if no survObj is provided, it is assumed that all individuals survive to the time of breeding
    if (is.null(survObj)) survObj <- makeSurvObj(data.frame(size=runif(21),surv=rep(1,21)),Formula=surv~size)
    # if no growObj is provided, it is assumed that all individuals retain the same size until the time of breeding		
    if (is.null(growObj)) growObj <- makeGrowthObj(data.frame(size=seq(21),sizeNext=seq(21)),Formula=sizeNext~size,Family="poisson")
    #print ("Warning: in the current version of IPMpack, makeIPMFmatrix still ignores the growObj you provided for your post-breeding F matrix. This will be included in a later version. Survival until breeding is already included in this version.")
    ##NOTE that the condition is necessary because you might ONLY have discrete offspring
    # in which case correction makes no sense
    if (fecObj@offspringSplitter$continuous>0) { 
      tmp <- t(outer(X=y,Y=y,.fecPostCensusInteger,
              cov=chosenCov,fecObj=fecObj, growObj=growObj,
              survObj=survObj, offspringObj= offspringObj))
    }
    
    
  }
  
  get.matrix <- to.cont <- tmp
  
  #discrete classes
  nDisc <- length(fecObj@offspringSplitter)-1
  namesDiscrete <- "NA"
  if (nDisc>0) {
    namesDiscrete <- colnames(fecObj@offspringSplitter[1:nDisc])
    to.discrete <- matrix(0,nrow=nDisc,ncol=nBigMatrix)
    for (i in 1:nDisc) {
      if (length(which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1))>1) {
        to.discrete[i,] <- apply(.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[2]][which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1),],2,prod)*unlist(fecObj@offspringSplitter[namesDiscrete[i]])
      } else {
        to.discrete[i,] <- .fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[2]][which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1),]*unlist(fecObj@offspringSplitter[namesDiscrete[i]])
      }
    }
    from.discrete <- matrix(0,ncol=nDisc,nrow=nDisc+nBigMatrix)
    if (names(fecObj@fecByDiscrete)[1]!="NA.") {
      if (sum(names(fecObj@fecByDiscrete)!=namesDiscrete)>0) stop ("Error - the names of the discrete classes as you provided for the data.frame fecByDiscrete are not 100% the same discrete class names in your data.frame offspringSplitter. They should also be in alphabetical order.")
      for (j in 1:nDisc) {
        for (i in 1:nDisc) {
          from.discrete[i,j] <- unlist(fecObj@offspringSplitter[namesDiscrete[i]]*fecObj@fecByDiscrete[namesDiscrete[j]])
        }
      }
      if (sum(fecObj@fecByDiscrete)>0&fecObj@offspringSplitter["continuous"]>0) {
        print ("WARNING - number and sizes of offspring produced by individuals in discrete classes cannot be calculated yet. The Fmatrix contains zeros instead. Only solution at this point: change the F matrix yourself afterwards.")
      }
    }
    get.matrix <- cbind(from.discrete,rbind(to.discrete,to.cont))
  }
  
  #print(colSums(get.matrix))
  
  #warning about negative numbers
  if (min(get.matrix)<0) { 
    print("Warning: fertility values < 0 exist in matrix, consider transforms. Negative values set to zero") 
    get.matrix[get.matrix<0] <- 0
  }
  
  rc <- new("IPMmatrix",
      nDiscrete = nDisc,
      nEnvClass = 1, 
      nBigMatrix = nBigMatrix,
      nrow = 1*nBigMatrix+nDisc,
      ncol =1*nBigMatrix+nDisc,
      meshpoints = y,
      env.index = rep(1:nEnvClass,each=nBigMatrix),
      names.discrete=namesDiscrete)
  rc[,] <-get.matrix   
  
  return(rc)
}

# =============================================================================
# =============================================================================
#Function creates a combo of F.IPMs for a chosen
#size range, with range env categories and transitions 
#between them
#and growth and survival objects
#
#Parameters -
#   nEnvClass - the number of env classes, defaults to 2, should match dim envMatrix
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   envMatrix - a matrix describing transiions between env
#   fecObj - a fecundity object
#   integrateType - NOT YET IMPLEMENTED
#
#Returns -
#  an IPM object


makeCompoundFmatrix <- function(nEnvClass = 2,
    nBigMatrix = 50,
    minSize = -1,
    maxSize = 50,
    envMatrix,
    fecObj,
    integrateType="midpoint",
    correction="none",
    preCensus=TRUE,
    survObj=NULL,
    growObj=NULL, offspringObj=NULL) {
  
  #warnings...
  if (nEnvClass!=nrow(envMatrix)) {
    print(paste("Dim of envMatrix not equal to nEnvClass. Adjusted to",nrow(envMatrix)))
    nEnvClass <- nrow(envMatrix)
  }
  
  # boundary points b and mesh points y
  b<-minSize+c(0:nBigMatrix)*(maxSize-minSize)/nBigMatrix;
  y<-0.5*(b[1:nBigMatrix]+b[2:(nBigMatrix+1)]);
  
  # step size for mid point rule, see equations 4 and 5
  h<-y[2]-y[1]
  
  #establish how how many discrete classes there are
  if (ncol(fecObj@offspringSplitter)>1) nDisc <- ncol(fecObj@offspringSplitter)-1 else nDisc <- 0
  
  #indexes for slotting in IPMs
  indexes <- rep(1:nEnvClass,each=(nBigMatrix+nDisc))
  #print(indexes)
  
  
  #megamatrix
  megamatrix <- matrix(0,(nBigMatrix+nDisc)*nEnvClass,(nBigMatrix+nDisc)*nEnvClass) 
  
  #loop over habitats / environments
  for (k in 1:nEnvClass) {
    
    get.matrix <- makeIPMFmatrix(nEnvClass =1,
        nBigMatrix = nBigMatrix,
        minSize = minSize,
        maxSize = maxSize,
        chosenCov = data.frame(covariate=as.factor(k)),
        fecObj=fecObj,
        integrateType=integrateType,
        correction=correction,preCensus=preCensus,
        survObj=survObj,
        growObj=growObj, offspringObj= offspringObj)
    
    #print(range(get.matrix))
    #print(get.matrix[1:5,1:5])
    
    # transit them
    subset <- c(1:nEnvClass)[envMatrix[,k]>0]                 
    for (j in subset) {
      megamatrix[indexes==j,indexes==k] <- get.matrix@.Data*envMatrix[j,k]; 
    }
    
  }
  
  #warning about negative numbers should appear from makeIPMFmatrix
  
  
  rc <- new("IPMmatrix",
      nEnvClass = nEnvClass, 
      nBigMatrix = nBigMatrix,
      nrow = nEnvClass*(nBigMatrix+nDisc),
      ncol =nEnvClass*(nBigMatrix+nDisc),
      meshpoints = y,
      env.index = rep(1:nEnvClass,each=nBigMatrix),
      names.discrete=get.matrix@names.discrete)
  
  
  rc[,] <- megamatrix
  
  return(rc) 
}



# =============================================================================
# =============================================================================
#Function creates a single F.IPM (fecundity transitions only)
#for a chosen size range, env category (single one at a time)
#and survival and fecundity objects (assume survival 
#precedes growth; could do otherwise...)
#
#Parameters -
#   nEnvClass - the number of env classes, cannot be !=1, defaults to 1 
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   chosenCov - current level of covariate
#   fecObj - a fecundity object
#   integrateType - options include "midpoint" "cumul" 
#   etc...
#Returns -
#  an IPM object


makeIPMCmatrix <- function(clonalObj,
    nEnvClass = 1,
    nBigMatrix = 50,
    minSize = -1,
    maxSize = 50,
    chosenCov = data.frame(covariate=1),
    integrateType="midpoint",
    correction="none", 
    preCensus=TRUE,
    survObj=NULL,
    growObj=NULL, offspringObj= NULL) {
  
  rc <- makeIPMFmatrix(fecObj=clonalObj,
      nEnvClass=nEnvClass,
      nBigMatrix=nBigMatrix,
      minSize=minSize,
      maxSize=maxSize,
      chosenCov=chosenCov,
      integrateType=integrateType,
      correction=correction, 
      preCensus=preCensus,
      survObj=survObj,
      growObj=growObj, offspringObj= offspringObj)
  
  return(rc)
}

# =============================================================================
# =============================================================================
#Function creates a combo of F.IPMs for a chosen
#size range, with range env categories and transitions 
#between them
#and growth and survival objects
#
#Parameters -
#   nEnvClass - the number of env classes, defaults to 2, should match dim envMatrix
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   envMatrix - a matrix describing transiions between env
#   fecObj - a fecundity object
#   integrateType - NOT YET IMPLEMENTED
#
#Returns -
#  an IPM object


makeCompoundCmatrix <- function(nEnvClass = 2,
    nBigMatrix = 50,
    minSize = -1,
    maxSize = 50,
    envMatrix,
    clonalObj,
    integrateType="midpoint",
    correction="none",
    preCensus=TRUE,
    survObj=NULL,
    growObj=NULL, offspringObj= NULL) {
  
  rc <- makeCompoundFmatrix(nEnvClass = nEnvClass,
      nBigMatrix = nBigMatrix,
      minSize = minSize,
      maxSize = maxSize,
      envMatrix = envMatrix,
      fecObj = clonalObj,
      integrateType=integrateType,
      correction=correction,
      preCensus=preCensus,
      survObj=survObj,
      growObj=growObj, offspringObj= offspringObj)
  
  return(rc) 
}


# =============================================================================
# =============================================================================
#Function creates a combo of T.IPMs for a chosen
#size range, with range env categories and transitions 
#between them
#and growth and survival objects
#
#Parameters -
#   nEnvClass - the number of env classes, defaults to 2, should match dim envMatrix
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   envMatrix - a matrix describing transiions between env
#   growObj - a growth object
#   survObj - a survival object
#   discreteTrans - discrete transition object, defaults to 1
#   integrateType - what type of integration?
#                    current options are "midpoint" (using pdf)
#                    or "cumul" (using cdf)
#   correction - do you not want to correct for integration errors ('none')
#                or correct by multiplying each column by a constant ('constant')
#
#Returns -
#  an IPM object


makeCompoundPmatrix <- function(nEnvClass = 2,
    nBigMatrix = 50,
    minSize = -1,
    maxSize = 50,
    envMatrix,
    growObj,
    survObj,
    discreteTrans=1,
    integrateType="midpoint",
    correction="none") {
  
  
  #warnings...
  if (nEnvClass!=nrow(envMatrix)) {
    print(paste("Dim of envMatrix not equal to nEnvClass. Adjusted to",nrow(envMatrix)))
    nEnvClass <- nrow(envMatrix)
  }
  
  # boundary points b and mesh points y
  b<-minSize+c(0:nBigMatrix)*(maxSize-minSize)/nBigMatrix;
  y<-0.5*(b[1:nBigMatrix]+b[2:(nBigMatrix+1)]);
  
  # step size for mid point rule, see equations 4 and 5
  h<-y[2]-y[1]
  
  #establish how how many discrete classes there are	
  if (class(discreteTrans) == "discreteTrans") nDisc <- ncol(discreteTrans@meanToCont) else nDisc <- 0
  
  
  #indexes for slotting in IPMs
  indexes <- rep(1:nEnvClass,each=(nBigMatrix+nDisc))
  
  #megamatrix
  megamatrix <- matrix(0,(nBigMatrix+nDisc)*nEnvClass,(nBigMatrix+nDisc)*nEnvClass) 
  
  #print(indexes)
  
  #loop over habitats / environments
  for (k in 1:nEnvClass) { 
    #IPM for individuals starting in env k
    get.matrix <- makeIPMPmatrix(nEnvClass = nEnvClass, 
        nBigMatrix = nBigMatrix, minSize = minSize, maxSize = maxSize, 
        chosenCov = data.frame(covariate = as.factor(k)), growObj=growObj, survObj=survObj, 
        discreteTrans = discreteTrans, 
        integrateType = integrateType, correction = correction)		
    
    # transit them
    subset <- c(1:nEnvClass)[envMatrix[,k]>0]
    for (j in subset) { 
      megamatrix[indexes==j,indexes==k] <- get.matrix*envMatrix[j,k]; 
    }
    
  }
  
  nmes<-""
  
  rc <- new("IPMmatrix",
      nEnvClass = nEnvClass, 
      nBigMatrix = nBigMatrix,
      nrow = nEnvClass*(nBigMatrix + nDisc),
      ncol = nEnvClass*(nBigMatrix + nDisc),
      meshpoints = y,
      env.index = rep(1:nEnvClass, each = nBigMatrix,
          names.discrete = nmes))
  
  rc[,] <- megamatrix
  
  return(rc) 
}

# =============================================================================
# =============================================================================
#Class for the matrix that holds the IPM
#setClass("DiscreteMatrix",
#		representation(nDiscrete = "numeric", #number of discrete classes
#				nEnvClass = "numeric", #number of covariate levels
#				nBigMatrix = "numeric", #the resolution of the IPM
#				meshpoints = "numeric",
#				env.index = "numeric",
#				names.discrete = "character"),
#		contains="matrix")



#Function creates a single T.IPM (survival and growth
#transitions only) for a chosen size range, env category
#(single one at a time) and growth and survival objects
#
#Parameters -
#   nEnvClass - the number of env classes, cannot be !=1, defaults to 1 
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   chosenCov - current level of covariate - can be vector for continuous
#                 temporal stochasticity case
#   growObj - a growth object
#   survObj - a survival object
#   discreteTrans - discrete transition object, defaults to 1
#   integrateType - what type of integration?
#                    current options are "midpoint" (using pdf)
#                    or "cumul" (using cdf)
#   correction - do you not want to correct for integration errors ('none')
#                or correct by multiplying each column by a constant ('constant')

#
#Returns -
#  an IPM object (with or without discrete classes)
makeIPMPmatrix <- function (nEnvClass = 1, nBigMatrix = 50, minSize = -1, maxSize = 50, 
    chosenCov = data.frame(covariate = 1), growObj, survObj, 
    discreteTrans = 1, integrateType = "midpoint", correction = "none") {
  
  if (class(growObj) == "growthObjPois" | class(growObj) =="growthObjNegBin") 
    print("warning: IPMs not appropriate with discrete growth processes")
  b <- minSize + c(0:nBigMatrix) * (maxSize - minSize)/nBigMatrix
  y <- 0.5 * (b[1:nBigMatrix] + b[2:(nBigMatrix + 1)])
  h <- y[2] - y[1]
  
  
  if (integrateType == "midpoint") {
    get.matrix <- t(outer(y, y, growSurv, cov = chosenCov, 
            growthObj = growObj, survObj = survObj)) * h
  }
  
  if (integrateType == "cumul") {
    get.matrix.cum <- t(outer(y, b, growthCum, cov = chosenCov, 
            growthObj = growObj))
    get.matrix <- get.matrix.cum[2:(nBigMatrix + 1), ] - 
        get.matrix.cum[1:nBigMatrix, ]
    #remove because alden...
    #get.matrix[nBigMatrix, nBigMatrix] <- get.matrix[nBigMatrix, 
    #			nBigMatrix] + (1 - sum(get.matrix[, nBigMatrix]))
    get.matrix <- t(t(get.matrix) * surv(size = y, cov = chosenCov, 
            survObj = survObj))
  }
  
  if (correction == "constant") {
    nvals <- colSums(get.matrix,na.rm=TRUE)
    loc0 <- which(nvals == 0 , arr.ind = TRUE)
    if (length(loc0) > 0) {
      print("warnings - columns that sum to 0 or that have NAs - assuming survival is along the diagonal; plot your Pmatrix to check it")
      get.matrix[,loc0] <- 0
      get.matrix[cbind(loc0, loc0)] <- surv(size = y[loc0], 
          cov = chosenCov, survObj = survObj)
    }
    nvals <- colSums(get.matrix,na.rm=TRUE)
    get.matrix <- t((t(get.matrix)/nvals) * surv(size = y, 
            cov = chosenCov, survObj = survObj))
  }
  
  if (correction=="discretizeExtremes"){
    tooLow <- growthCum(y,b[1], cov = chosenCov, 
        growthObj = growObj)
    tooHigh <- 1-growthCum(y,b[length(b)],cov = chosenCov, 
        growthObj = growObj)
    get.matrix[1,] <- get.matrix[1,]+tooLow * surv(size = y, cov = chosenCov, survObj = survObj)
    get.matrix[nBigMatrix,] <- get.matrix[nBigMatrix,]+tooHigh * surv(size = y, cov = chosenCov, survObj = survObj)
  }
  
  
  rc <- new("IPMmatrix", nDiscrete = 0, nEnvClass = 1, nBigMatrix = nBigMatrix, 
      nrow = 1 * nBigMatrix, ncol = 1 * nBigMatrix, meshpoints = y, 
      env.index = rep(1:nEnvClass, each = nBigMatrix), names.discrete = "")
  rc[, ] <- get.matrix
  
  if (class(discreteTrans) == "discreteTrans") {
    nDisc <- ncol(discreteTrans@meanToCont)
    moveToDiscrete <- predict(discreteTrans@moveToDiscrete, 
        data.frame(size = y, size2 = (y * y)), type = "response")
    cont.to.cont <- get.matrix * matrix(1 - moveToDiscrete, 
        nrow = nBigMatrix, ncol = nBigMatrix, byrow = TRUE)
    disc.to.disc <- discreteTrans@discreteTrans[1:nDisc, 1:nDisc]
    disc.to.cont <- matrix(0, ncol = nDisc, nrow = nBigMatrix)
    cont.to.disc <- matrix(0, nrow = nDisc, ncol = nBigMatrix)
    for (j in 1:nDisc) {
      tmp <- dnorm(y, discreteTrans@meanToCont[j], discreteTrans@sdToCont[j]) * h
      if (correction == "constant") tmp <- tmp/sum(tmp)
      tmp[which(is.na(tmp))] <- 0
      disc.to.cont[, j]  <- discreteTrans@discreteTrans["continuous", j] * tmp
      if (discreteTrans@discreteTrans[j,"continuous"]>0) {
        cont.to.disc[j, ] <- surv(y, chosenCov, survObj) * moveToDiscrete * 
            discreteTrans@discreteTrans[j,"continuous"] / sum(discreteTrans@discreteTrans[1:nDisc,"continuous"]) 
      }
    }
    
    
    get.disc.matrix <- rbind(cbind(disc.to.disc, cont.to.disc), 
        cbind(disc.to.cont, cont.to.cont))
    rc <- new("IPMmatrix", nDiscrete = nDisc, nEnvClass = 1, 
        nBigMatrix = nBigMatrix, nrow = 1 * nBigMatrix + 
            nDisc, ncol = 1 * nBigMatrix + nDisc, meshpoints = y, 
        env.index = rep(1:nEnvClass, each = nBigMatrix), 
        names.discrete = rownames(discreteTrans@discreteTrans)[1:nDisc])
    rc[, ] <- get.disc.matrix
  }
  return(rc)
}

# =============================================================================
# =============================================================================
##same function... but for discrete (no integration!)  ####
###
### NOTE ASSUMES TRANSITIONS OUT ARE DISCRETE - THIS MEANS THAT SDTOCONT
### IS THE NEGBINON ISSUE    #########
##
makeIntegerPmatrix <- function (nEnvClass = 1, 
    meshpoints=1:20,
    chosenCov = data.frame(covariate = 1), 
    growObj, survObj, 
    discreteTrans = 1) {
  
  y <- meshpoints
  nBigMatrix <- length(y)
  
  get.matrix <- t(outer(y, y, growSurv, cov = chosenCov, 
          growthObj = growObj, survObj = survObj))
  
  
  rc <- new("IPMmatrix", nDiscrete = 0, nEnvClass = 1, nBigMatrix = nBigMatrix, 
      nrow = 1 * nBigMatrix, ncol = 1 * nBigMatrix, meshpoints = y, 
      env.index = rep(1:nEnvClass, each = nBigMatrix), names.discrete = "")
  rc[, ] <- get.matrix
  
  if (class(discreteTrans) == "discreteTrans") {
    nDisc <- ncol(discreteTrans@meanToCont)
    moveToDiscrete <- predict(discreteTrans@moveToDiscrete, 
        data.frame(size = y, size2 = (y * y)), type = "response")
    cont.to.cont <- get.matrix * matrix(1 - moveToDiscrete, 
        nrow = nBigMatrix, ncol = nBigMatrix, byrow = TRUE)
    disc.to.disc <- discreteTrans@discreteTrans[1:nDisc, 1:nDisc]
    disc.to.cont <- matrix(0, ncol = nDisc, nrow = nBigMatrix)
    cont.to.disc <- matrix(0, nrow = nDisc, ncol = nBigMatrix)
    for (j in 1:nDisc) {
      if (discreteTrans@distToCont=="poisson") 
        tmp <- dpois(y, discreteTrans@meanToCont[j]) 
      if (discreteTrans@distToCont=="negBin") 
        tmp <- dnbinom(y, mu=discreteTrans@meanToCont[j], size=discreteTrans@thetaToCont[j]) 
      tmp[which(is.na(tmp))] <- 0
      disc.to.cont[, j] <- discreteTrans@discreteTrans["continuous", j] * tmp
      if (sum(discreteTrans@discreteTrans[j,"continuous"]>0)) {
        cont.to.disc[j, ] <- surv(y, chosenCov, survObj) * moveToDiscrete * 
            discreteTrans@discreteTrans[j,"continuous"] / sum(discreteTrans@discreteTrans[1:nDisc,"continuous"]) 
      }
    }	
    get.disc.matrix <- rbind(cbind(disc.to.disc, cont.to.disc), 
        cbind(disc.to.cont, cont.to.cont))
    
    rc <- new("DiscreteMatrix", nDiscrete = nDisc, nEnvClass = 1, 
        nBigMatrix = nBigMatrix, nrow = 1 * nBigMatrix + 
            nDisc, ncol = 1 * nBigMatrix + nDisc, meshpoints = y, 
        env.index = rep(1:nEnvClass, each = nBigMatrix), 
        names.discrete = rownames(discreteTrans@discreteTrans)[1:nDisc])
    rc[, ] <- get.disc.matrix
  }
  return(rc)
}

# =============================================================================
# =============================================================================
#Function creates a single F.IPM (fecundity transitions only)
#for a chosen size range, env category (single one at a time)
#and survival and fecundity objects (assume survival 
#precedes growth; could do otherwise...)
#
#Parameters -
#   nEnvClass - the number of env classes, cannot be !=1, defaults to 1 
#   nBigMatrix - the number of size bins in the model
#   minSize - lower end of desired size range
#   maxSize - upper end of desired size range
#   chosenCov - current level of covariate
#   fecObj - a fecundity object
#   integrateType - options include "midpoint" "cumul" 
#   etc...
#Returns -
#  an IPM object

makeIPMFmatrix <- function(fecObj,
    nEnvClass = 1,
    nBigMatrix = 50,
    minSize = -1,
    maxSize = 50,
    chosenCov = data.frame(covariate=1),
    integrateType="midpoint",
    correction="none",
    preCensus=TRUE,
    survObj=NULL,
    growObj=NULL,
    offspringObj=NULL) {
  
  # boundary points b and mesh points y
  b<-minSize+c(0:nBigMatrix)*(maxSize-minSize)/nBigMatrix;
  y<-0.5*(b[1:nBigMatrix]+b[2:(nBigMatrix+1)]);
  
  # step size for mid point rule, see equations 4 and 5
  h<-y[2]-y[1]
  
  fecObj@fecConstants[is.na(fecObj@fecConstants)] <- 1
  
  tmp <- matrix(0,ncol=length(y),nrow=length(y))
  
  # 1. pre-census
  if (preCensus) { 
    #print("here")
    ##NOTE that the condition is necessary because you might ONLY have discrete offspring
    # in which case correction makes no sense
    if (integrateType=="midpoint"&fecObj@offspringSplitter$continuous>0) { 
      tmp <- t(outer(X=y,Y=y,.fecPreCensus,cov=chosenCov,fecObj=fecObj, offspringObj=offspringObj))*h 
    }
    ##NOTE that the condition is necessary because you might ONLY have discrete offspring
    # in which case correction makes no sense
    if (integrateType=="cumul"&fecObj@offspringSplitter$continuous>0) {
      #offspring extremes (pnorm) 
      tmp.cum <- t(outer(X=y,Y=b,.offspringCum,cov=chosenCov,
              fecObj=fecObj,offspringObj=offspringObj))
      tmp <- tmp.cum[2:(nBigMatrix+1),]-tmp.cum[1:nBigMatrix,]
      #put in seed production
      tmp <- t(t(tmp)*.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]])      
    }
    
    ##NOTE that the condition is necessary because you might ONLY have discrete offspring
    # in which case correction makes no sense
    if (correction=="constant"&fecObj@offspringSplitter$continuous>0) {
      # in this case, column sums should equal raw fecundity
      correction.here <-.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]]/colSums(tmp)
      tmp <- t(t(tmp)*correction.here)
    }
    ##NOTE that the condition is necessary because you might ONLY have discrete offspring
    # in which case correction makes no sense
    if (correction=="discretizeExtremes"&fecObj@offspringSplitter$continuous>0){
      tooLow <- .offspringCum(x=y,y=b[1], cov = chosenCov, 
          fecObj = fecObj,offspringObj=offspringObj)
      tooHigh <- 1-.offspringCum(x=y,y=b[length(b)],cov = chosenCov, 
          fecObj = fecObj,offspringObj=offspringObj)
      tmp[1,] <- tmp[1,]+tooLow
      tmp[nBigMatrix,] <- tmp[nBigMatrix,]+tooHigh
    }
    
    # 2. post-census
  } else {
    # if no survObj is provided, it is assumed that all individuals survive to the time of breeding
    if (is.null(survObj)) survObj <- makeSurvObj(data.frame(size=runif(21),surv=rep(1,21)),Formula=surv~size)
    # if no growObj is provided, it is assumed that all individuals retain the same size until the time of breeding		
    if (is.null(growObj)) growObj <- makeGrowthObj(data.frame(size=seq(21),sizeNext=seq(21)),Formula=sizeNext~size)
    #print ("Warning: in the current version of IPMpack, makeIPMFmatrix still ignores the growObj you provided for your post-breeding F matrix. This will be included in a later version. Survival until breeding is already included in this version.")
    ##NOTE that the condition is necessary because you might ONLY have discrete offspring
    # in which case correction makes no sense
    if (integrateType=="midpoint"&fecObj@offspringSplitter$continuous>0) { 
      tmp <- t(outer(X=y,Y=y,.fecPostCensus,
              cov=chosenCov,fecObj=fecObj, growObj=growObj,
              survObj=survObj,offspringObj=offspringObj))*h 
    }
    ##NOTE that the condition is necessary because you might ONLY have discrete offspring
    # in which case correction makes no sense
    if (integrateType=="cumul"&fecObj@offspringSplitter$continuous>0) {
      #make the extreme bins offspring matrix
      tmp.cum <- t(outer(X=y,Y=b,.offspringCum,cov=chosenCov,
              fecObj=fecObj,offspringObj=offspringObj))
      tmpBabies <- tmp.cum[2:(nBigMatrix+1),]-tmp.cum[1:nBigMatrix,]
      
      #make the extreme bins growth matrix
      tmp.cum <- t(outer(X=y,Y=b,growthCum,cov=chosenCov,
              growObj=growObj,offspringObj=offspringObj))
      tmpGrowth <- tmp.cum[2:(nBigMatrix+1),]-tmp.cum[1:nBigMatrix,]
      tmpGrowth[nBigMatrix,nBigMatrix] <- tmpGrowth[nBigMatrix,nBigMatrix]+
          (1-sum(tmpGrowth[,nBigMatrix]))
      
      #put in survival and seed production
      tmp <- t(t(tmpBabies*tmpGrowth)*surv(size=y,cov=chosenCov,survObj=survObj)*.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]])
      
    }
    
    ##NOTE that the condition is necessary because you might ONLY have discrete offspring
    # in which case correction makes no sense
    if (correction=="constant"&fecObj@offspringSplitter$continuous>0) {
      # in this case, column sums should equal raw fecundity * survival
      correction.here <-.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[1]]*surv(size=y,cov=chosenCov,survObj=survObj)/colSums(tmp)
      tmp <- t(t(tmp)*correction.here)
    }
    ##NOTE that the condition is necessary because you might ONLY have discrete offspring
    # in which case correction makes no sense
    if (correction=="discretizeExtremes"&fecObj@offspringSplitter$continuous>0){
      tooLow <- .offspringCum(x=y,y=b[1], cov = chosenCov, 
          fecObj = fecObj,offspringObj=offspringObj)
      tooHigh <- 1-.offspringCum(x=y,y=b[length(b)],cov = chosenCov, 
          fecObj = fecObj,offspringObj=offspringObj)
      tmp[1,] <- tmp[1,]+tooLow
      tmp[nBigMatrix,] <- tmp[nBigMatrix,]+tooHigh
    }
  }
  
  get.matrix <- to.cont <- tmp
  
  #discrete classes
  nDisc <- length(fecObj@offspringSplitter)-1
  namesDiscrete <- "NA"
  if (nDisc>0) {
    namesDiscrete <- colnames(fecObj@offspringSplitter[1:nDisc])
    to.discrete <- matrix(0,nrow=nDisc,ncol=nBigMatrix)
    for (i in 1:nDisc) {
      if (length(which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1))>1) {
        to.discrete[i,] <- apply(.fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[2]][which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1),],2,prod)*unlist(fecObj@offspringSplitter[namesDiscrete[i]])
      } else {
        to.discrete[i,] <- .fecRaw(x=y,cov=chosenCov,fecObj=fecObj)[[2]][which(fecObj@vitalRatesPerOffspringType[,namesDiscrete[i]]==1),]*unlist(fecObj@offspringSplitter[namesDiscrete[i]])
      }
    }
    from.discrete <- matrix(0,ncol=nDisc,nrow=nDisc+nBigMatrix)
    if (names(fecObj@fecByDiscrete)[1]!="NA.") {
      if (sum(names(fecObj@fecByDiscrete)!=namesDiscrete)>0) stop ("Error - the names of the discrete classes as you provided for the data.frame fecByDiscrete are not 100% the same discrete class names in your data.frame offspringSplitter. They should also be in alphabetical order.")
      for (j in 1:nDisc) {
        for (i in 1:nDisc) {
          from.discrete[i,j] <- unlist(fecObj@offspringSplitter[namesDiscrete[i]]*fecObj@fecByDiscrete[namesDiscrete[j]])
        }
      }
      if (sum(fecObj@fecByDiscrete)>0&fecObj@offspringSplitter["continuous"]>0) {
        print ("WARNING - number and sizes of offspring produced by individuals in discrete classes cannot be calculated yet. The Fmatrix contains zeros instead. Only solution at this point: change the F matrix yourself afterwards.")
      }
    }
    get.matrix <- cbind(from.discrete,rbind(to.discrete,to.cont))
  }
  
  #print(colSums(get.matrix))
  
  #warning about negative numbers
  if (min(get.matrix)<0) { 
    print("Warning: fertility values < 0 exist in matrix, consider transforms. Negative values set to zero") 
    get.matrix[get.matrix<0] <- 0
  }
  
  rc <- new("IPMmatrix",
      nDiscrete = nDisc,
      nEnvClass = 1, 
      nBigMatrix = nBigMatrix,
      nrow = 1*nBigMatrix+nDisc,
      ncol =1*nBigMatrix+nDisc,
      meshpoints = y,
      env.index = rep(1:nEnvClass,each=nBigMatrix),
      names.discrete=namesDiscrete)
  rc[,] <-get.matrix   
  
  return(rc)
}


# ===========================================================================
# ===========================================================================
# function to make an IPM from P, F, and C matrices that has all the same slot definitions as the P matrix
makeIPMmatrix = function(Pmatrix,Fmatrix,Cmatrix=NULL){
  if(is.null(Cmatrix)){ #make an empty placeholder
    Cmatrix=Fmatrix
    Cmatrix@.Data=0*Fmatrix@.Data
  }
  IPMmatrix=Pmatrix
  IPMmatrix@.Data <- 	IPMmatrix@.Data + Fmatrix@.Data + Cmatrix@.Data
  return(IPMmatrix)
}

# ===========================================================================
# ===========================================================================

