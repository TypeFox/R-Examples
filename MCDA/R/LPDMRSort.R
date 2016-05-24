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

LPDMRSort <- function(performanceTable, categoriesLowerProfiles, criteriaWeights, criteriaMinMax, majorityThreshold, criteriaVetos = NULL, criteriaDictators = NULL, majorityRule = "", alternativesIDs = NULL, criteriaIDs = NULL, categoriesIDs = NULL){
  
  ## check the input data
  
  if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
    stop("wrong performanceTable, should be a matrix or a data frame")
  
  if (!(is.matrix(categoriesLowerProfiles)))
    stop("categoriesLowerProfiles should be a matrix")
  
  if (!(is.vector(criteriaMinMax)))
    stop("criteriaMinMax should be a vector")
  
  if (!(is.vector(criteriaWeights)))
    stop("criteriaWeights should be a vector")
  
  if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
    stop("alternativesIDs should be a vector")
  
  if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
    stop("criteriaIDs should be a vector")
  
  if (!(is.null(categoriesIDs) || is.vector(categoriesIDs)))
    stop("categoriesIDs should be a vector")
  
  if (!(is.null(criteriaVetos) || is.matrix(criteriaVetos)))
    stop("criteriaVetos should be a matrix")
  
  if (!(is.null(criteriaDictators) || is.matrix(criteriaDictators)))
    stop("criteriaDictators should be a matrix")
  
  if (!is.character(majorityRule))
    stop("majorityRule should be a string")
  else if (!(majorityRule %in% c("","V","D","v","d","dV","Dv","dv")))
    stop("majorityRule needs to take values in {'','V','D','v','d','dV','Dv','dv'}")
  
  if (majorityRule %in% c("V","v","dV","Dv","dv") && is.null(criteriaVetos))
    stop("majorityRule requires non-NULL criteriaVetos")
  
  if (majorityRule %in% c("D","d","dV","Dv","dv") && is.null(criteriaDictators))
    stop("majorityRule requires non-NULL criteriaDictators")
  
  ## filter the data according to the given alternatives and criteria
  
  if (!is.null(alternativesIDs)){
    performanceTable <- performanceTable[alternativesIDs,]
  } 
  
  if (!is.null(criteriaIDs)){
    performanceTable <- performanceTable[,criteriaIDs]
    criteriaWeights <- criteriaWeights[criteriaIDs]
    criteriaMinMax <- criteriaMinMax[criteriaIDs]
    categoriesLowerProfiles <- categoriesLowerProfiles[,criteriaIDs]
  }
  
  if ((!is.null(criteriaIDs)) && (!is.null(criteriaVetos))){
    criteriaVetos <- criteriaVetos[,criteriaIDs]  
  }
  
  if ((!is.null(criteriaIDs)) && (!is.null(criteriaDictators))){
    criteriaDictators <- criteriaDictators[,criteriaIDs]  
  }
  
  if ((!is.null(categoriesIDs)) && (!is.null(criteriaVetos))){
    criteriaVetos <- criteriaVetos[categoriesIDs,]
  }
    
  if (!is.null(categoriesIDs)){
    categoriesLowerProfiles <- categoriesLowerProfiles[categoriesIDs,]
  }
  
  if ((!is.null(categoriesIDs)) && (!is.null(criteriaDictators))){
    criteriaDictators <- criteriaDictators[categoriesIDs,]
  }  
  
  # data is filtered, check for some data consistency
  
  # if there are less than 2 criteria or 2 alternatives, there is no MCDA problem
  
  if (is.null(dim(performanceTable))) 
    stop("less than 2 criteria or 2 alternatives")
  
  # -------------------------------------------------------
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  numCat <- dim(categoriesLowerProfiles)[1]
  
  # -------------------------------------------------------
  
  outranking <- function(alternativePerformances, profilePerformances, criteriaWeights, criteriaMinMax, majorityThreshold, profileCriteriaVetos=NULL, profileCriteriaDictators=NULL, majorityRule = ""){
    localConcordance <- rep(0,numCrit)
    veto <- 0
    dictator <- 0
    for (i in 1:numCrit)
    {
      if (criteriaMinMax[i] == "min")
      {
        if (alternativePerformances[i] %<=% profilePerformances[i])
          localConcordance[i] = 1
        if (!is.null(profileCriteriaVetos))
        {
          if (!is.na(profileCriteriaVetos[i]))
          {
            if (alternativePerformances[i] %>=% profileCriteriaVetos[i])
              veto = 1 
          }
        }
        if (!is.null(profileCriteriaDictators))
        {
          if (!is.na(profileCriteriaDictators[i]))
          {
            if (alternativePerformances[i] %<=% profileCriteriaDictators[i])
              dictator = 1 
          }
        }
      }
      else
      {
        if (alternativePerformances[i] %>=% profilePerformances[i])
          localConcordance[i] = 1
        if (!is.null(profileCriteriaVetos))
        {
          if (!is.na(profileCriteriaVetos[i]))
          {
            if (alternativePerformances[i] %<=% profileCriteriaVetos[i])
              veto = 1 
          }
        }
        if (!is.null(profileCriteriaDictators))
        {
          if (!is.na(profileCriteriaDictators[i]))
          {
            if (alternativePerformances[i] %>=% profileCriteriaDictators[i])
              dictator = 1 
          }
        }
      }
    }
    
    concordance = sum(localConcordance*criteriaWeights)
    
    if(majorityRule == "")
    {
      if(!(concordance %>=% majorityThreshold))
        return(FALSE)
      else
        return(TRUE)
    }
    else if(majorityRule == "V")
    {
      if ((veto == 1) || !(concordance %>=% majorityThreshold))
        return(FALSE)
      else 
        return(TRUE)
    }
    else if(majorityRule == "D")
    {
      if ((dictator == 1) || (concordance %>=% majorityThreshold))
        return(TRUE)
      else 
        return(FALSE)
    }
    else if(majorityRule == "v")
    {
      if ((veto == 1 && dictator == 0) || !(concordance %>=% majorityThreshold))
        return(FALSE)
      else
        return(TRUE)
    }
    else if(majorityRule == "d")
    {
      if ((dictator == 1 && veto == 0) || (concordance %>=% majorityThreshold))
        return(TRUE)
      else 
        return(FALSE)
    }
    else if(majorityRule == "Dv")
    {
      if ((dictator == 1) || (concordance %>=% majorityThreshold && veto == 0))
        return(TRUE)
      else 
        return(FALSE)
    }
    else if(majorityRule == "dV")
    {
      if ((veto == 1) || (!(concordance %>=% majorityThreshold) && dictator == 0))
        return(FALSE)
      else 
        return(TRUE)
    }
    else if(majorityRule == "dv")
    {
      if ((dictator == 1 && veto == 0) || (concordance %>=% majorityThreshold && veto == 0) || (concordance %>=% majorityThreshold && veto == 1 && dictator == 1))
        return(TRUE)
      else 
        return(FALSE)
    }
  }
  
  assignments <- c()
  
  for (i in 1:numAlt)
  {
    categoryNotFound <- TRUE
    k <- 1
    while ((categoryNotFound) && k<=numCat-1)
    {
      profileCriteriaVetos <- NULL
      if (!is.null(criteriaVetos))
        profileCriteriaVetos <- criteriaVetos[k,]
      
      profileCriteriaDictators <- NULL
      if (!is.null(criteriaDictators))
        profileCriteriaDictators <- criteriaDictators[k,]

      if(outranking(performanceTable[i,],categoriesLowerProfiles[k,], criteriaWeights, criteriaMinMax, majorityThreshold, profileCriteriaVetos = profileCriteriaVetos, profileCriteriaDictators = profileCriteriaDictators, majorityRule = majorityRule))
      {
        category <- k
        categoryNotFound <- FALSE
      }
      else
        k<-k+1
    }
    assignments <- c(assignments, rownames(categoriesLowerProfiles)[k])
  }
  
  names(assignments) <- rownames(performanceTable)
  
  return(assignments)
}
