#############################################################################
#
# Copyright Patrick Meyer and Sébastien Bigaret, 2013
#
# Contributors:
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

applyPiecewiseLinearValueFunctionsOnPerformanceTable <- function(valueFunctions, performanceTable, alternativesIDs = NULL, criteriaIDs = NULL){
  
  ## check the input data
  
  if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
    stop("wrong performanceTable, should be a matrix or a data frame")
  
  if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
    stop("alternatives IDs should be in a vector")
  
  if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
    stop("criteria IDs should be in a vector")
  
  ## filter the data according to the given criteria and alternatives
  
  if (!is.null(criteriaIDs)){
    valueFunctions <- valueFunctions[criteriaIDs]
    performanceTable <- performanceTable[,criteriaIDs]
  }
  
  if (!is.null(alternativesIDs)){
    performanceTable <- performanceTable[alternativesIDs,]
  } 
  
  # -------------------------------------------------------
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  # -------------------------------------------------------
  
  ## check if performanceTable values lie in the ranges of the valueFunctions values
  
  test <- TRUE
  i<-1
  
  while(test && (i<=numCrit)){
    mini <- min(performanceTable[,i])
    maxi <- max(performanceTable[,i])
    
    vFmini <- min(valueFunctions[[i]]["x",])
    vFmaxi <- max(valueFunctions[[i]]["x",])
    
    if (!((mini >= vFmini) && (maxi <= vFmaxi)))
      test<-FALSE
    
    i<-i+1
  }
  
  if (!test)
    stop("performanceTable ranges do not lie in valueFunctions ranges")
  
  
  ## reorder points of value functions to make them increasing w.r.t. y
  
  for (i in 1:numCrit){
    valueFunctions[[i]] <- valueFunctions[[i]][,order(valueFunctions[[i]][2,])]
  }
  
  ## determine if criteria are to be minimized or maximized according to the value functions
  
  criteriaMinMax<-c()
  
  for (i in 1:numCrit){
    if (all(diff(valueFunctions[[i]][1,]) %>=% 0)){
      criteriaMinMax <- c(criteriaMinMax, "max")
    }
    else if (all(diff(valueFunctions[[i]][1,]) %<=% 0)){
      criteriaMinMax <- c(criteriaMinMax, "min")
    }
    else stop("non monotonic value function")
  }
  
  # determine first how many breakpoints in each value function
  
  criteriaNumberOfBreakPoints <- c()
  
  for (i in 1:numCrit)
    criteriaNumberOfBreakPoints <- c(criteriaNumberOfBreakPoints, dim(valueFunctions[[i]])[2])
  
  criteriaBreakPoints <- list()
  
  for (i in 1:numCrit){
    tmp <- valueFunctions[[i]]["x",]
    criteriaBreakPoints <- c(criteriaBreakPoints,list(tmp))
  }
  
  names(criteriaBreakPoints) <- colnames(performanceTable)
  
  # -------------------------------------------------------
  
  # a is a matrix decomposing the alternatives in the break point space
  
  a<-matrix(0,nrow=numAlt, ncol=(sum(criteriaNumberOfBreakPoints)))
  
  for (n in 1:numAlt){
    for (m in 1:numCrit){
      if (length(which(performanceTable[n,m]==criteriaBreakPoints[[m]]))!=0){
        # then we have a performance value which is on a breakpoint
        j<-which(performanceTable[n,m]==criteriaBreakPoints[[m]])
        if (m==1)
          pos <- j
        else
          pos<-sum(criteriaNumberOfBreakPoints[1:(m-1)])+j
        a[n,pos] <- 1
      }
      else{
        # then we have value which needs to be approximated by a linear interpolation
        # let us first search the lower and upper bounds of the interval of breakpoints around the value
        if (criteriaMinMax[m] == "min"){
          j<-which(performanceTable[n,m]>criteriaBreakPoints[[m]])[1]-1
        }
        else{
          j<-which(performanceTable[n,m]<criteriaBreakPoints[[m]])[1]-1			
        }
        if (m==1)
          pos <- j
        else
          pos<-sum(criteriaNumberOfBreakPoints[1:(m-1)])+j
        a[n,pos] <- 1-(performanceTable[n,m]-criteriaBreakPoints[[m]][j])/(criteriaBreakPoints[[m]][j+1] - criteriaBreakPoints[[m]][j])
        a[n,pos+1] <- (performanceTable[n,m]-criteriaBreakPoints[[m]][j])/(criteriaBreakPoints[[m]][j+1] - criteriaBreakPoints[[m]][j])
      }
    }
  }
  
  # -------------------------------------------------------
  
  normalizedPerformanceTable <- matrix(nrow=numAlt, ncol=numCrit)
  
  for (m in 1:numCrit){
    if (m==1){
      start = 1
      end = criteriaNumberOfBreakPoints[m]
    }
    else
    {
      start = sum(criteriaNumberOfBreakPoints[1:(m-1)]) + 1
      end = sum(criteriaNumberOfBreakPoints[1:(m-1)]) + criteriaNumberOfBreakPoints[m]
    }
    normalizedPerformanceTable[,m] <- a[,start:end]%*%valueFunctions[[m]]["y",]
  }
  rownames(normalizedPerformanceTable) <- rownames(performanceTable)
  
  colnames(normalizedPerformanceTable) <- colnames(performanceTable)
  
  return(normalizedPerformanceTable)
}
