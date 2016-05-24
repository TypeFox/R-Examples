#############################################################################
#
# Copyright Alexandru Olteanu, Patrick Meyer, and Sébastien Bigaret, 2015
#
# Contributors:
#   Alexandru Olteanu <al.olteanu@gmail.com>
#   Patrick Meyer <patrick.meyer@telecom-bretagne.eu>
#   Sébastien Bigaret <sebastien.bigaret@telecom-bretagne.eu>
#  	
# This software, MCDA, is a package for the R statistical software which 
# allows to use MCDA algorithms and methods. 
# 
# This software is governed by the CeCILL license (v2) under French law
# and abiding by the rules of distribution of free software. You can
# use, modify and/ or redistribute the software under the terms of the
# CeCILL license as circulated by CEA, CNRS and INRIA at the following
# URL "http://www.cecill.info".
# 
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
#		
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#		
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##############################################################################

LPDMRSortInferenceExact <- function(performanceTable, assignments, categoriesRanks, criteriaMinMax, majorityRule = "", readableWeights = FALSE, readableProfiles = FALSE, minmaxLPD = FALSE, alternativesIDs = NULL, criteriaIDs = NULL){
  
  ## check the input data
  if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
    stop("wrong performanceTable, should be a matrix or a data frame")
  
  if (!(is.vector(assignments)))
    stop("assignments should be a vector")
  
  if (!(is.vector(categoriesRanks)))
    stop("categoriesRanks should be a vector")
  
  if (!(is.vector(criteriaMinMax)))
    stop("criteriaMinMax should be a vector")
  
  if (!is.character(majorityRule))
    stop("majorityRule should be a boolean")
  else if (!(majorityRule %in% c("","V","D","v","d","dV","Dv","dv")))
    stop("majorityRule needs to take values in {'','V','D','v','d','dV','Dv','dv'}")
  
  if (!is.logical(readableWeights))
    stop("readableWeights should be a boolean")
  
  if (!is.logical(readableProfiles))
    stop("readableProfiles should be a boolean")
  
  if (!is.logical(minmaxLPD))
    stop("minmaxLPD should be a boolean")
  
  if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
    stop("alternativesIDs should be a vector")
  
  if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
    stop("criteriaIDs should be a vector")
  
  ## filter the data according to the given alternatives and criteria
  
  if (!is.null(alternativesIDs)){
    performanceTable <- performanceTable[alternativesIDs,]
    assignments <- assignments[alternativesIDs]
  } 
  
  if (!is.null(criteriaIDs)){
    performanceTable <- performanceTable[,criteriaIDs]
    criteriaMinMax <- criteriaMinMax[criteriaIDs]
  }
  
  # data is filtered, check for some data consistency
  
  # if there are less than 2 criteria or 2 alternatives, there is no MCDA problem
  
  if (is.null(dim(performanceTable))) 
    stop("less than 2 criteria or 2 alternatives")
  
  # -------------------------------------------------------
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  numCat <- length(categoriesRanks)
  
  tempPath <- tempdir()
  
  # get model file depending on function options

  modelfilename <- paste("MRSort", c("","V","D","DV1","DV2","DV3","DV4","DV5")[match(majorityRule,c("","V","D","v","d","dV","Dv","dv"))], "InferenceModel", sep = "")
  
  if(readableWeights || readableProfiles)
  {
    modelfilename <- paste(modelfilename, "Spread", sep = "")
    if(readableWeights)
      modelfilename <- paste(modelfilename, "Weights", sep = "")
    if(readableProfiles)
      modelfilename <- paste(modelfilename, "Profiles", sep = "")
  }
  if(minmaxLPD)
    modelfilename <- paste(modelfilename, "LPD", sep = "")
  
  modelfilename <- paste(modelfilename, ".gmpl", sep = "")
  
  modelFile <- system.file("extdata", modelfilename, package="MCDA")
  
  dataFile <- tempfile()
  
  file.copy(modelFile, dataFile)
  
  sink(dataFile, append=TRUE)
  
  cat("data;\n")
  
  cat("param X := ")
  cat(numAlt)
  cat(";\n\n")
  
  cat("param F := ")
  cat(numCrit)
  cat(";\n\n")
  
  cat("param Fdir := \n")
  for (i in 1:numCrit){
    cat(i)
    cat("\t")
    if (criteriaMinMax[i]=="min")
      cat("-1")
    else
      cat("1")
    if (i!=numCrit)
      cat("\n")
    else
      cat(";\n\n")
  }
  
  cat("param Fmin :=\n")
  for (i in 1:numCrit){
    cat(i)
    cat("\t")
    cat(apply(performanceTable, 2, min)[i])
    if (i!=numCrit)
      cat("\n")
    else
      cat(";\n\n") 
  }
  
  cat("param Fmax :=\n")
  for (i in 1:numCrit){
    cat(i)
    cat("\t")
    cat(apply(performanceTable, 2, max)[i])
    if (i!=numCrit)
      cat("\n")
    else
      cat(";\n\n") 
  }
  
  cat("param K := ")
  cat(numCat)
  cat(";\n\n")
  
  cat("param A:=\n")
  for (i in 1:numAlt){
    cat(i)
    cat("\t")
    cat(categoriesRanks[assignments[i]])
    if (i!=numAlt)
      cat("\n")
    else
      cat(";\n\n") 
  }
  
  cat("param PTx : ")
  cat(1:numCrit)
  cat(" := \n")
  for (i in 1:numAlt){
    cat(i)
    cat("\t")
    cat(performanceTable[i,])
    if (i!=numAlt)
      cat("\n")
    else
      cat(";\n\n")
  }
  
  cat("param gamma:=0.001;\n")
  
  cat("end;\n")
  sink()
  
  lp<-initProbGLPK()
  
  tran<-mplAllocWkspGLPK()
  
  setMIPParmGLPK(PRESOLVE, GLP_ON)
  
  termOutGLPK(GLP_OFF)
  
  out<-mplReadModelGLPK(tran, dataFile, skip=0)
  
  if (is.null(out))
    out <- mplGenerateGLPK(tran)
  else 
    stop(return_codeGLPK(out))
  
  if (is.null(out))
    mplBuildProbGLPK(tran,lp)
  else 
    stop(return_codeGLPK(out))
  
  solveMIPGLPK(lp)
  
  if(mipStatusGLPK(lp)==5){
    
    mplPostsolveGLPK(tran, lp, sol = GLP_MIP)
    
    solution <- mipColsValGLPK(lp)
    
    varnames <- c()
    
    for (i in 1:length(solution))
      varnames <- c(varnames,getColNameGLPK(lp,i))
    
    lambda <- solution[varnames=="lambda"]
    
    weightsnames <- c()
    
    for (i in 1:numCrit)
    {
      weightsnames <- c(weightsnames,paste("w[",i,"]",sep=""))
    }
    
    weights <- c()
    
    for (i in 1:numCrit)
      weights <- c(weights,solution[varnames==weightsnames[i]])
    
    names(weights) <- colnames(performanceTable)
    
    ptknames <- matrix(nrow=numCat,ncol=numCrit)
    
    for (i in 2:(numCat+1)){
      for (j in 1:numCrit)
      {
        ptknames[i-1,j] <- paste("PTk[",i,",",j,"]",sep="")
      }
    }
    
    profilesPerformances <- matrix(nrow=numCat,ncol=numCrit)
    
    for (i in 1:numCat){
      for (j in 1:numCrit)
        profilesPerformances[i,j] <- solution[varnames==ptknames[i,j]]
    }
    
    rownames(profilesPerformances) <- names(categoriesRanks)
    colnames(profilesPerformances) <- colnames(performanceTable)
    
    vetoPerformances <- NULL
    
    if(majorityRule %in% c("V","v","d","dV","Dv","dv"))
    {
      ptvnames <- matrix(nrow=numCat,ncol=numCrit)
      
      for (i in 2:(numCat+1)){
        for (j in 1:numCrit)
        {
          ptvnames[i-1,j] <- paste("PTv[",i,",",j,"]",sep="")
        }
      }
      
      vetoPerformances <- matrix(nrow=numCat,ncol=numCrit)
      
      for (i in 1:numCat){
        for (j in 1:numCrit)
          vetoPerformances[i,j] <- solution[varnames==ptvnames[i,j]]
      }
      
      rownames(vetoPerformances) <- names(categoriesRanks)
      colnames(vetoPerformances) <- colnames(performanceTable)
    }
    
    dictatorPerformances <- NULL
    
    if(majorityRule %in% c("D","v","d","dV","Dv","dv"))
    {
      ptdnames <- matrix(nrow=numCat,ncol=numCrit)
      
      for (i in 2:(numCat+1)){
        for (j in 1:numCrit)
        {
          ptdnames[i-1,j] <- paste("PTd[",i,",",j,"]",sep="")
        }
      }
      
      dictatorPerformances <- matrix(nrow=numCat,ncol=numCrit)
      
      for (i in 1:numCat){
        for (j in 1:numCrit)
          dictatorPerformances[i,j] <- solution[varnames==ptdnames[i,j]]
      }
      
      rownames(dictatorPerformances) <- names(categoriesRanks)
      colnames(dictatorPerformances) <- colnames(performanceTable)
    }
    
    return(list(lambda = lambda, weights = weights, profilesPerformances = profilesPerformances, vetoPerformances = vetoPerformances, dictatorPerformances = dictatorPerformances))
    
  }
  else
    return(NULL)
}
