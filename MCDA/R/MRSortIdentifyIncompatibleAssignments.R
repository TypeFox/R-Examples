#############################################################################
#
# Copyright Alexandru Olteanu, Patrick Meyer and Sébastien Bigaret, 2015
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

MRSortIdentifyIncompatibleAssignments <- function(performanceTable, assignments, categoriesRanks, criteriaMinMax, veto = FALSE, incompatibleSetsLimit = 100, largerIncompatibleSetsMargin = 0, alternativesIDs = NULL, criteriaIDs = NULL){
  
  ## check the input data
  if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
    stop("wrong performanceTable, should be a matrix or a data frame")
  
  if (!(is.vector(assignments)))
    stop("assignments should be a vector")
  
  if (!(is.vector(categoriesRanks)))
    stop("categoriesRanks should be a vector")
  
  if (!(is.vector(criteriaMinMax)))
    stop("criteriaMinMax should be a vector")
  
  if (!is.logical(veto))
    stop("veto should be a boolean")
  
  if (!is.numeric(incompatibleSetsLimit))
    stop("incompatibleSetsLimit should be numeric")
  else if (incompatibleSetsLimit%%1!=0)
    stop("incompatibleSetsLimit should be an integer")
  else if (incompatibleSetsLimit<=0)
    stop("incompatibleSetsLimit should be strictly pozitive")
  
  if (!is.numeric(largerIncompatibleSetsMargin))
    stop("largerIncompatibleSetsMargin should be numeric")
  else if (largerIncompatibleSetsMargin%%1!=0)
    stop("largerIncompatibleSetsMargin should be an integer")
  else if (largerIncompatibleSetsMargin<0)
    stop("largerIncompatibleSetsMargin should be pozitive")
  
  if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
    stop("alternativesIDs should be a vector")
  
  if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
    stop("criteriaIDs should be a vector")
  
  ## filter the data according to the given alternatives and criteria
  
  if (!is.null(alternativesIDs)){
    performanceTable <- performanceTable[alternativesIDs,]
    assignments <- assignments[alternativesIDs]
  }
  else
    alternativesIDs = rownames(performanceTable)
  
  if (!is.null(criteriaIDs)){
    performanceTable <- performanceTable[,criteriaIDs]
    criteriaMinMax <- criteriaMinMax[criteriaIDs]
  }
  else
    criteriaIDs = colnames(performanceTable)
  
  # data is filtered, check for some data consistency
  
  # if there are less than 2 criteria or 2 alternatives, there is no MCDA problem
  
  if (is.null(dim(performanceTable))) 
    stop("less than 2 criteria or 2 alternatives")
  
  # -------------------------------------------------------
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  numCat <- length(categoriesRanks)
  
  tempPath <- tempdir()
  
  # get data content that remains the same for all following linear program executions
  
  datacontent <- paste("data;\nparam X := ", numAlt, ";\n\nparam F := ", numCrit, ";\n\nparam Fdir := \n", sep = "")
  
  for (i in 1:numCrit){
    datacontent <- paste(datacontent, i, "\t", sep = "")
    if (criteriaMinMax[i]=="min")
      datacontent <- paste(datacontent, "-1", sep = "")
    else
      datacontent <- paste(datacontent, "1", sep = "")
    if (i!=numCrit)
      datacontent <- paste(datacontent, "\n", sep = "")
    else
      datacontent <- paste(datacontent, ";\n\n", sep = "")
  }
  
  datacontent <- paste(datacontent, "param Fmin :=\n", sep = "")
  
  for (i in 1:numCrit){
    datacontent <- paste(datacontent, i, "\t", apply(performanceTable, 2, min)[i], sep = "")
    if (i!=numCrit)
      datacontent <- paste(datacontent, "\n", sep = "")
    else
      datacontent <- paste(datacontent, ";\n\n", sep = "")
  }
  
  datacontent <- paste(datacontent, "param Fmax :=\n", sep = "")
  
  for (i in 1:numCrit){
    datacontent <- paste(datacontent, i, "\t", apply(performanceTable, 2, max)[i], sep = "")
    if (i!=numCrit)
      datacontent <- paste(datacontent, "\n", sep = "")
    else
      datacontent <- paste(datacontent, ";\n\n", sep = "")
  }
  
  datacontent <- paste(datacontent, "param K :=", numCat, ";\n\n", sep = "")
  
  datacontent <- paste(datacontent, "param A:=\n", sep = "")
  
  for (i in 1:numAlt){
    datacontent <- paste(datacontent, i, "\t", categoriesRanks[assignments[i]], sep = "")
    if (i!=numAlt)
      datacontent <- paste(datacontent, "\n", sep = "")
    else
      datacontent <- paste(datacontent, ";\n\n", sep = "")
  }
  
  datacontent <- paste(datacontent, "param PTx : ", sep = "")
  for(i in 1:numCrit)
    datacontent <- paste(datacontent, i, sep = " ")
  datacontent <- paste(datacontent, ":= \n", sep = "")
  
  
  for (i in 1:numAlt){
    datacontent <- paste(datacontent, i, "\t", sep = "")
    for (j in 1:numCrit){
      datacontent <- paste(datacontent, performanceTable[i,j], sep = "")
      if (j!=numCrit)
        datacontent <- paste(datacontent, " ", sep = "")
    }
    if (i!=numAlt)
      datacontent <- paste(datacontent, "\n", sep = "")
    else
      datacontent <- paste(datacontent, ";\n\n", sep = "")
  }
  
  datacontent <- paste(datacontent, "param gamma:=0.0001;\n", sep = "")
  
  # get first model file
  
  modelFile <- system.file("extdata","MRSortIdentifyMinimalInvalidAssignmentsSet.gmpl", package="MCDA")
  if(veto)
    modelFile <- system.file("extdata","MRSortVIdentifyMinimalInvalidAssignmentsSet.gmpl", package="MCDA")
  
  # write data file
  
  dataFile <- tempfile()
  
  file.copy(modelFile, dataFile)
  
  sink(dataFile, append=TRUE)
  
  cat(datacontent)
  
  cat("end;\n")
  
  sink()
  
  # init and run linear program
  
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
    
    # get variable names
    
    varnames <- c()
    
    for (i in 1:length(solution))
      varnames <- c(varnames,getColNameGLPK(lp,i))
    
    # get size of minimal incompatible assignments set and one such set
    
    minIncompatibleSetsSize <- 0
    
    incompatibleSet <- c()
    
    for (i in 1:numAlt)
    {
      if(solution[varnames==paste("OnOff[",i,"]",sep="")] == 1)
      {
        incompatibleSet <- c(incompatibleSet,alternativesIDs[i])
        minIncompatibleSetsSize <- minIncompatibleSetsSize + 1
      }
    }
    
    incompatibleSets <- list(incompatibleSet)
    
    # if there are no incompatible sets return the empty set
    
    if(minIncompatibleSetsSize == 0)
      return(incompatibleSets)
    
    # get second model file
    
    modelFile <- system.file("extdata","MRSortIdentifyInvalidAssignmentsSet.gmpl", package="MCDA")
    if(veto)
      modelFile <- system.file("extdata","MRSortVIdentifyInvalidAssignmentsSet.gmpl", package="MCDA")
    
    # create new data content
    
    datacontent2a <- "param PrevOnOff : "
    for(i in 1:numAlt)
      datacontent2a <- paste(datacontent2a, i, sep = " ")
    datacontent2a <- paste(datacontent2a, ":= \n1\t", sep = "")
    for(i in 1:numAlt)
      datacontent2a <- paste(datacontent2a, solution[varnames==paste("OnOff[",i,"]",sep="")], sep = " ")
    
    datacontent2b <- paste("param PrevOnOffLimit := \n1\t ", minIncompatibleSetsSize, sep ="")

    # iterate through acceptes sizes for incompatible assignment sets
    
    incompatibleSetSize <- minIncompatibleSetsSize
    
    while(incompatibleSetSize <= minIncompatibleSetsSize + largerIncompatibleSetsMargin)
    {
      # break if we've retrieved the desired number of incompatible sets
      
      if(length(incompatibleSets) >= incompatibleSetsLimit)
        break
      
      repeat{
        # write data file
        
        dataFile <- tempfile()
        
        file.copy(modelFile, dataFile)
        
        sink(dataFile, append=TRUE)
        
        cat(datacontent)
        
        cat("param invalid:=")
        cat(incompatibleSetSize)
        cat(";\n")
        
        cat("param Y:=")
        cat(length(incompatibleSets))
        cat(";\n")
        
        cat(datacontent2a)
        cat(";\n\n")
        
        cat(datacontent2b)
        cat(";\n\n")
        
        cat("end;\n")
        
        sink()
        
        # init and run linear program
        
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
        
        if(mipStatusGLPK(lp)==5)
        {
          # found new incompatible assignments set
          
          mplPostsolveGLPK(tran, lp, sol = GLP_MIP)
          
          solution <- mipColsValGLPK(lp)
          
          varnames <- c()
          
          for (i in 1:length(solution))
            varnames <- c(varnames,getColNameGLPK(lp,i))
          
          # get incompatible assignments set
          
          incompatibleSet <- c()
          
          for (i in 1:numAlt)
            if(solution[varnames==paste("OnOff[",i,"]",sep="")] == 1)
              incompatibleSet <- c(incompatibleSet,alternativesIDs[i])
          
          # add set
          
          incompatibleSets <- c(incompatibleSets, list(incompatibleSet))
          
          # update data content
          
          datacontent2a <- paste(datacontent2a, "\n", length(incompatibleSets), "\t", sep = "")
          for(i in 1:numAlt)
            datacontent2a <- paste(datacontent2a, solution[varnames==paste("OnOff[",i,"]",sep="")], sep = " ")
          
          datacontent2b <- paste(datacontent2b, "\n", length(incompatibleSets), "\t", incompatibleSetSize, sep ="")
          
        }
        else
          break
      }
      
      # increase size of incompatible sets
      
      incompatibleSetSize <- incompatibleSetSize + 1
    }
    
    return(incompatibleSets)
    
  }
  else
    return(NULL)
}
