#############################################################################
#
#  This file is a part of the R package "RoughSets".
#
#  Author: Lala Septem Riza and Andrzej Janusz
#  Supervisors: Chris Cornelis, Francisco Herrera, Dominik Slezak and Jose Manuel Benitez
#  Copyright (c):
#       DiCITS Lab, Sci2s group, DECSAI, University of Granada and
#       Institute of Mathematics, University of Warsaw
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################

# an auxiliary function for discretizing a single attribute
# @param vec a vector of data
# @param cuts a vector of cut values
applyDiscretization <- function(vec, cuts, isNominal) {
   if(!isNominal) vec = cut(vec, c(-Inf,cuts,Inf),
                            right=TRUE, include.lowest = TRUE,
                            ordered_result = TRUE)
   return(vec)
}

# a function for computing cuts of the "quantile-based" discretization of a single attribute into a given number of intervals
# @param vec a vector of data
# @param n a number
discretize.quantiles <- function(vec, n) {
  uniqueValues = unique(vec)
  if (length(uniqueValues) < n) {
		if(length(uniqueValues) == 1) cutVec = NULL
		else cutVec = uniqueValues
	}
	else {
		cutVec = stats::quantile(vec, (1:(n-1))/n)
		cutVec = unique(cutVec)
	}

  return(cutVec)
}

# a function for computing cuts of the "equal interval size" discretization of a single attribute into a given number of intervals
# @param vec a vector of data
# @param n a number
discretize.equal.intervals <- function(vec, n) {
  attrRange = range(vec)
  if(diff(attrRange) > 0 && n > 1) {
		cutVec = seq(attrRange[1], attrRange[2], length.out = n + 1)
  }
  else {
		cutVec = NULL
  }

  return(cutVec[-c(1,length(cutVec))])
}

# a function for computing cuts using the maximum discernibility heuristic (the global approach)
global.discernibility <- function(vecList, cutCandidatesVecList, decVec, nOfCuts,
                                 nAttrs = length(vecList), minIntSupport = 0, ...) {
  attrCount = length(vecList)
  cutVecList = list()
  rmVecList = list()
  notRemovedVecIndicator = sapply(cutCandidatesVecList, function(x) return(length(x) > 0))
  cutVecList[1:attrCount] = list(numeric())

  INDclasses = list(1:length(decVec))
  minIntervalSize = ceiling(length(decVec)*minIntSupport)
  scrHistVec = rep(0, attrCount)
#  maxDiscernPairs = conflictsCouner(decVec)
  i = 0
  numOfChosenCuts = 0
  nDecisions = length(unique(decVec))
  endFlag = FALSE

  while(numOfChosenCuts < nOfCuts && !endFlag) {
    i = i + 1
    rmVecList[1:attrCount] = list(integer())
    candidateVecIdx = which(notRemovedVecIndicator)
    if(nAttrs < length(candidateVecIdx)) {
      attrSampleIdx = sample(candidateVecIdx, nAttrs, replace = FALSE)
    } else {
      attrSampleIdx = candidateVecIdx
    }
    tmpINDclassesVec = unlist(INDclasses)
    tmpObjIdxLengths = sapply(INDclasses, length)

    bestCutsList = mapply(evaluateCuts, vecList[attrSampleIdx], cutCandidatesVecList[attrSampleIdx],
                          MoreArgs = list(decVec = decVec[tmpINDclassesVec],
                                          nOfDec = nDecisions,
                                          INDclassesVec = tmpINDclassesVec,
                                          INDclassesSizes = tmpObjIdxLengths,
                                          minIntervalSize = minIntervalSize),
                          SIMPLIFY = F)

    maxCutScoreVec = sapply(bestCutsList, function(x) x$maxTPtoFP)
    maxCutIdxVec = sapply(bestCutsList, function(x) x$maxIdx)
    rmVecList[attrSampleIdx] = lapply(bestCutsList, function(x) x$rmVec)
    rm(tmpINDclassesVec, tmpObjIdxLengths, bestCutsList)

    maxScr = max(maxCutScoreVec)
    tmpIdxVec = which(maxCutScoreVec == maxScr)
    bestAttrIdx = tmpIdxVec[which.max(scrHistVec[attrSampleIdx[tmpIdxVec]])]
    numOfChosenCuts = numOfChosenCuts + length(bestAttrIdx)
    scrHistVec[attrSampleIdx] = maxCutScoreVec
    chosenCutIdx = maxCutIdxVec[bestAttrIdx]
    rm(maxCutScoreVec, maxCutIdxVec, tmpIdxVec)

    if(maxScr == 0 || numOfChosenCuts >= nOfCuts) {
      endFlag = TRUE
      if(numOfChosenCuts >= nOfCuts)
        cutVecList[[attrSampleIdx[bestAttrIdx]]] = c(cutVecList[[attrSampleIdx[bestAttrIdx]]],
                                                     cutCandidatesVecList[[attrSampleIdx[bestAttrIdx]]][chosenCutIdx])
    }
    else  {
      cutVecList[[attrSampleIdx[bestAttrIdx]]] = c(cutVecList[[attrSampleIdx[bestAttrIdx]]],
                                                   cutCandidatesVecList[[attrSampleIdx[bestAttrIdx]]][chosenCutIdx])
      rmVecList[[attrSampleIdx[bestAttrIdx]]] = c(rmVecList[[attrSampleIdx[bestAttrIdx]]],
                                                  chosenCutIdx)

      tmpDiscretizedVec = as.integer(vecList[[attrSampleIdx[bestAttrIdx]]] >=
                                       cutCandidatesVecList[[attrSampleIdx[bestAttrIdx]]][chosenCutIdx])
      newINDclasses = unlist(lapply(INDclasses, splitINDclass, tmpDiscretizedVec), recursive = FALSE)
      tmpClassesToRmIdx = which(sapply(newINDclasses, function(x, decV) length(unique(decV[x])) == 1,
                                       as.character(decVec)))
      if(length(tmpClassesToRmIdx) > 0) newINDclasses = newINDclasses[-tmpClassesToRmIdx]

      INDclasses = newINDclasses
      if(length(INDclasses) > 0)  maxDiscernPairs = sum(sapply(INDclasses, function(x, decV) conflictsCouner(decV[x]), decVec))
      else maxDiscernPairs = 0
      rm(newINDclasses, tmpClassesToRmIdx, tmpDiscretizedVec)

      for(j in attrSampleIdx) {
        if(length(rmVecList[[j]]) > 0)  {
          cutCandidatesVecList[[j]] = cutCandidatesVecList[[j]][-rmVecList[[j]]]
          if(length(cutCandidatesVecList[[j]]) == 0) notRemovedVecIndicator[j] = FALSE
        }
      }

      if(maxDiscernPairs == 0 || max(sapply(INDclasses,length)) < 2*minIntervalSize)  endFlag = TRUE
    }
  }

  return(cutVecList)
}


# a function for computing cuts using the local discernibility heuristic (the local approach)
local.discernibility <- function(vec, cutCandidatesVec, decVec,
                                 nDecisions, nOfCuts = 2, minIntSupport = 0) {

  cutVec = numeric()

  INDclasses = list(1:length(decVec))
  minIntervalSize = ceiling(length(decVec)*minIntSupport)
  numOfChosenCuts = 0
  endFlag = FALSE

  while(!endFlag) {
    tmpINDclassesVec = unlist(INDclasses)
    tmpObjIdxLengths = sapply(INDclasses, length)

    bestCut = evaluateCuts(vec, cutCandidatesVec,
                           decVec = decVec[tmpINDclassesVec],
                           nOfDec = nDecisions,
                           INDclassesVec = tmpINDclassesVec,
                           INDclassesSizes = tmpObjIdxLengths,
                           minIntervalSize = minIntervalSize)

    maxScr = bestCut$maxTPtoFP
    chosenCutIdx = bestCut$maxIdx
    rm(tmpINDclassesVec, tmpObjIdxLengths, bestCut)

    numOfChosenCuts = numOfChosenCuts + 1

    if(maxScr == 0 || numOfChosenCuts >= nOfCuts) {
      endFlag = TRUE
      if(numOfChosenCuts >= nOfCuts)
        cutVec = c(cutVec, cutCandidatesVec[chosenCutIdx])
    } else  {
      cutVec = c(cutVec, cutCandidatesVec[chosenCutIdx])
      rmCutIdx = chosenCutIdx

      tmpDiscretizedVec = as.integer(vec >= cutCandidatesVec[chosenCutIdx])
      newINDclasses = unlist(lapply(INDclasses, splitINDclass, tmpDiscretizedVec), recursive = FALSE)

      tmpClassesToRmIdx = which(sapply(newINDclasses, function(x, decV) length(unique(decV[x])) == 1,
                                       as.character(decVec)))
      if(length(tmpClassesToRmIdx) > 0) newINDclasses = newINDclasses[-tmpClassesToRmIdx]

      INDclasses = newINDclasses
      if(length(INDclasses) > 0)  {
        maxDiscernPairs = sum(sapply(INDclasses, function(x, decV) conflictsCouner(decV[x]), decVec))
      } else maxDiscernPairs = 0
      rm(newINDclasses, tmpClassesToRmIdx, tmpDiscretizedVec)

      cutCandidatesVec = cutCandidatesVec[-rmCutIdx]

      if(maxDiscernPairs == 0 || max(sapply(INDclasses,length)) < 2*minIntervalSize)  endFlag = TRUE
    }
  }

  return(cutVec)
}


#  an auxiliary function for computation of a number of conflicts with a decision vector
conflictsCouner <- function(decisionVector)  {
   decisionDistrib = as.numeric(table(decisionVector))
   return(as.numeric(sum(as.numeric(sum(decisionDistrib) - decisionDistrib) * decisionDistrib)))
}

# an auxiliary function for evaluation of a set of candidate cuts using the number of conflicts
# it uses a C++ code to speed up the computations
evaluateCuts <- function(numVec, cutsCandidates, decVec, nOfDec,
                        INDclassesVec, INDclassesSizes, minIntervalSize) {

  bestCut = .C("chooseBestCutC", k = as.integer(length(cutsCandidates)),
               cutCandidates = as.double(cutsCandidates),
               N = as.integer(length(INDclassesVec)),
               vec = as.double(numVec),
               objectsIdx = as.integer(INDclassesVec),
               objectsIdxLengths = as.integer(INDclassesSizes),
               numOfInt = as.integer(length(INDclassesSizes)),
               decVec = as.integer(decVec),
               nOfDec = as.integer(nOfDec),
               attrType = as.integer(T),
               minIntervalSize = as.integer(minIntervalSize),
               rmVec = as.integer(rep(0, length(cutsCandidates))),
               idxVec = as.integer(0),
               maxTPtoFP = as.double(0.0))

  return(list(maxTPtoFP = as.numeric(bestCut$maxTPtoFP),
              maxIdx = as.integer(bestCut$idxVec),
              rmVec = which(bestCut$rmVec > 0)))
}

# an auxiliary function for construction of sets of candidate cuts using for a given attribute
# it uses a C++ code to speed up the computations
chooseCutCandidates <- function(attrVec, decVec)  {

  tmpIdx = order(attrVec)
  tmpCutCandidates = .C("chooseCutCandidatesC", vec = as.double(attrVec[tmpIdx]),
                        decVec = as.integer(decVec[tmpIdx]),
                        N = as.integer(length(decVec)),
                        candidatesIdx = as.integer(rep(FALSE,length(decVec)-1)),
                        candidates = as.double(rep(0.0,length(decVec)-1)))
  return(unique(tmpCutCandidates$candidates[which(as.logical(tmpCutCandidates$candidatesIdx))]))
}

