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

normalizePerformanceTable <- function(performanceTable, normalizationTypes=NULL, alternativesIDs = NULL, criteriaIDs = NULL){
  
  ## http://en.wikipedia.org/wiki/Feature_scaling
  
  ## check the input data
  
  if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
    stop("wrong performance table, should be a matrix or a data frame")
  
  if (!(is.null(normalizationTypes) || is.vector(normalizationTypes)))
    stop("normalizationTypes should be a vector")
  
  if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
    stop("alternatives IDs should be in a vector")
  
  if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
    stop("criteria IDs should be in a vector")
  
  ## filter the performance table and the criteria according to the given alternatives and criteria
  
  if (!is.null(alternativesIDs)) performanceTable <- performanceTable[alternativesIDs,]
  
  if (!is.null(criteriaIDs)) performanceTable <- performanceTable[,criteriaIDs]
  
  if (!is.null(criteriaIDs) && !is.null(normalizationTypes)) normalizationTypes <- normalizationTypes[criteriaIDs]
  
  for (i in 1:dim(performanceTable)[2]){
    if (normalizationTypes[i] == "percentageOfMax"){
      performanceTable[,i] <- percentageOfMax(performanceTable[,i])
    }
    else if (normalizationTypes[i] == "rescaling"){
      performanceTable[,i] <- rescaling(performanceTable[,i])
    }
    else if (normalizationTypes[i] == "standardization"){
      performanceTable[,i] <- standardization(performanceTable[,i])
    }
    else if (normalizationTypes[i] == "scaleToUnitLength"){
      performanceTable[,i] <- scaleToUnitLength(performanceTable[,i])
    }
  }

  return(performanceTable)
  
}

percentageOfMax <- function(data){
  max <- max(data)
  data <- data/max
  return(data)
}

rescaling <- function(data){
  max <- max(data)
  min <- min(data)
  data <- (data-min)/(max-min)
  return(data)
}

standardization <- function(data){
  mean <- mean(data)
  sd <- sd(data)
  data <- (data - mean)/sd
  return(data)
}

scaleToUnitLength <- function(data){
  norm <- sqrt(sum(data^2))
  data <- data/norm
  return(data)
}
