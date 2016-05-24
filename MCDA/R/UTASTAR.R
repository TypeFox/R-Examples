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

UTASTAR <- function(performanceTable, criteriaMinMax, criteriaNumberOfBreakPoints, epsilon, alternativesRanks = NULL, alternativesPreferences = NULL, alternativesIndifferences = NULL, criteriaLBs=NULL, criteriaUBs=NULL, alternativesIDs = NULL, criteriaIDs = NULL, kPostOptimality = NULL){
  
  ## check the input data
  
  if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
    stop("wrong performanceTable, should be a matrix or a data frame")
  
  if (!(is.null(alternativesRanks) || is.vector(alternativesRanks)))
    stop("alternativesRanks should be a vector")
  
  if (!(is.null(alternativesPreferences) || is.matrix(alternativesPreferences)))
    stop("alternativesPreferences should be a matrix")
  
  if (!(is.null(alternativesIndifferences) || is.matrix(alternativesIndifferences)))
    stop("alternativesIndifferences should be a matrix")
  
  if (is.null(alternativesRanks) && is.null(alternativesPreferences) && is.null(alternativesIndifferences))
    stop("at least one of alternativesRanks, alternativesPreferences or alternativesIndifferences should not be NULL")
  
  if (!is.null(alternativesRanks) && (!is.null(alternativesPreferences) | !is.null(alternativesIndifferences)))
    stop("alternativesRanks and one of alternativesPreferences or alternativesIndifferences cannot be simultaneously not NULL")
  
  if (!(is.vector(criteriaMinMax)))
    stop("criteriaMinMax should be a vector")
  
  if (!(is.vector(criteriaNumberOfBreakPoints)))
    stop("criteriaNumberOfBreakPoints should be a vector")
  
  if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
    stop("alternativesIDs should be in a vector")
  
  if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
    stop("criteriaIDs should be in a vector")
  
  if (!(is.null(criteriaLBs) || is.vector(criteriaLBs)))
    stop("criteriaLBs should be in a vector")
  
  if (!(is.null(criteriaUBs) || is.vector(criteriaUBs)))
    stop("criteriaUBs should be in a vector")
  
  ## filter the data according to the given alternatives and criteria
  ## in alternativesIDs and criteriaIDs
  
  if (!is.null(alternativesIDs)){
    performanceTable <- performanceTable[alternativesIDs,]
    if (!is.null(alternativesRanks))
      alternativesRanks <- alternativesRanks[alternativesIDs]
    if (!is.null(alternativesPreferences)){
      tmpIds <- intersect(alternativesPreferences, alternativesIDs)
      tmpMatrix <- c()
      for (i in 1:dim(alternativesPreferences)[1]){
        if (all(alternativesPreferences[i,] %in% tmpIds))
          tmpMatrix <- rbind(tmpMatrix,alternativesPreferences[i,])
      }
      alternativesPreferences <- tmpMatrix
    }
    if (!is.null(alternativesIndifferences)){
      tmpIds <- intersect(alternativesIndifferences, alternativesIDs)
      tmpMatrix <- c()
      for (i in 1:dim(alternativesIndifferences)[1]){
        if (all(alternativesIndifferences[i,] %in% tmpIds))
          tmpMatrix <- rbind(tmpMatrix,alternativesIndifferences[i,])
      }
      alternativesIndifferences <- tmpMatrix
    }
  } 
  
  if (!is.null(criteriaIDs)){
    criteriaMinMax <- criteriaMinMax[criteriaIDs]
    performanceTable <- performanceTable[,criteriaIDs]
    criteriaNumberOfBreakPoints <- criteriaNumberOfBreakPoints[criteriaIDs]
  }
  
  if (!is.null(criteriaIDs) && !is.null(criteriaUBs)){
    criteriaUBs <- criteriaUBs[criteriaIDs]
  }
  
  if (!is.null(criteriaIDs) && !is.null(criteriaLBs)){
    criteriaLBs <- criteriaLBs[criteriaIDs]
  }
  
  # only the alternatives which are in the ranking should be considered for the calculation
  # we therefore take the intersection between the alternatives present in the performance
  # table and those of the ranking
  
  if (!is.null(alternativesRanks)){
    reallyActiveAlternatives <- intersect(rownames(performanceTable),names(alternativesRanks))
    
    if (length(reallyActiveAlternatives) != 0){
      performanceTable <- performanceTable[reallyActiveAlternatives,]
      alternativesRanks <- alternativesRanks[reallyActiveAlternatives]
    } else {
      stop("alternatives of alternativesRanks are not compatible with those of performanceTable")
    }
  }
  
  if (!is.null(alternativesPreferences) || !is.null(alternativesIndifferences)){
    reallyActiveAlternatives <- intersect(rownames(performanceTable),rbind(alternativesPreferences,alternativesIndifferences))
    
    if (length(reallyActiveAlternatives) != 0){
      performanceTable <- performanceTable[reallyActiveAlternatives,]
      
      if (!is.null(alternativesPreferences)){
        tmpIds <- intersect(alternativesPreferences, reallyActiveAlternatives)
        tmpMatrix <- c()
        for (i in 1:dim(alternativesPreferences)[1]){
          if (all(alternativesPreferences[i,] %in% tmpIds))
            tmpMatrix <- rbind(tmpMatrix,alternativesPreferences[i,])
        }
        alternativesPreferences <- tmpMatrix
      }
      
      if (!is.null(alternativesIndifferences)){
        tmpIds <- intersect(alternativesIndifferences, reallyActiveAlternatives)
        tmpMatrix <- c()
        for (i in 1:dim(alternativesIndifferences)[1]){
          if (all(alternativesIndifferences[i,] %in% tmpIds))
            tmpMatrix <- rbind(tmpMatrix,alternativesIndifferences[i,])
        }
        alternativesIndifferences <- tmpMatrix
      }
      
    } else {
      stop("alternatives of alternativesPreferences or alternativesIndifferences are not compatible with those of performanceTable")
    }
    
  }
  
  # data is filtered, check for some data consistency
  
  # are the upper and lower bounds given in the function compatible with the data in the performance table ?
  if (!(is.null(criteriaUBs))){
    if (!all(apply(performanceTable,2,max)<=criteriaUBs))
      stop("performanceTable contains higher values than criteriaUBs")
  }
  
  if (!(is.null(criteriaLBs))){
    if (!all(apply(performanceTable,2,min)>=criteriaLBs))
      stop("performanceTable contains lower values than criteriaLBs")
  }
  
  if (!all(criteriaNumberOfBreakPoints >= 2))
    stop("in criteriaNumberOfBreakPoints there should at least be 2 breakpoints for each criterion")
  
  # if there are less than 2 criteria or 2 alternatives, there is no MCDA problem
  
  if (is.null(dim(performanceTable))) 
    stop("less than 2 criteria or 2 alternatives")
  
  
  # if there are no alternatives left in the ranking or the pairwise preferences
  # we stop here
  
  if (is.null(alternativesRanks) && is.null(alternativesPreferences) && is.null(alternativesIndifferences))
    stop("after filtering none of alternativesRanks, alternativesPreferences or alternativesIndifferences is not NULL")
  
  # -------------------------------------------------------
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  # -------------------------------------------------------
  
  criteriaBreakPoints <- list()
  
  for (i in 1:numCrit){
    
    tmp<-c()
    
    if (!is.null(criteriaLBs))
      mini <- criteriaLBs[i]
    else{
      mini <- min(performanceTable[,i])
    }
    
    if (!is.null(criteriaLBs))
      maxi <- criteriaUBs[i]
    else{
      maxi <- max(performanceTable[,i])
    }
    
    if (mini == maxi){
      # then there is only one value for that criterion, and the algorithm to build the linear interpolation
      # will not work correctly
      stop(paste("there is only one possible value left for criterion "),colnames(performanceTable)[i])
    }
    
    alphai <- criteriaNumberOfBreakPoints[i]
    
    options(digits=20)
    
    for (j in 1:alphai)
      tmp<-c(tmp,mini + (j-1)/(alphai-1) * (maxi - mini))
    
    # due to this formula, the minimum and maximum values might not be exactly the same than the real minimum and maximum values in the performance table
    # to be sure there is no rounding problem, we recopy these values in tmp (important for the later comparisons to these values)
    
    tmp[1] <- mini
    tmp[alphai] <- maxi
    
    # if the criterion has to be maximized, the worst value is in the first position
    # else, we sort the vector the other way around to have the worst value in the first position
    
    if (criteriaMinMax[i] == "min")
      tmp<-sort(tmp,decreasing=TRUE)
    criteriaBreakPoints <- c(criteriaBreakPoints,list(tmp))
  }
  
  names(criteriaBreakPoints) <- colnames(performanceTable)
  
  # -------------------------------------------------------
  
  # a is a matrix decomposing the alternatives in the break point space and adding the sigmaPlus and sigmaMinus columns
  
  a<-matrix(0,nrow=numAlt, ncol=(sum(criteriaNumberOfBreakPoints)+2*numAlt))
  
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
      # and now for sigmaPlus
      a[n,dim(a)[2]-2*numAlt+n] <- -1
      # and sigmaMinus
      a[n,dim(a)[2]-numAlt+n] <- +1
    }
  }
  
  # -------------------------------------------------------
  
  # the objective function : the first elements correspond to the ui's, the last one to the sigmas
  
  obj<-rep(0,sum(criteriaNumberOfBreakPoints))
  
  obj<-c(obj,rep(1,2*numAlt))
  
  # -------------------------------------------------------
  
  # we now build the part of the constraints matrix concerning the order / preferences / indifferences given by the decision maker
  
  preferenceConstraints<-matrix(nrow=0, ncol=sum(criteriaNumberOfBreakPoints)+2*numAlt)
  indifferenceConstraints <-matrix(nrow=0, ncol=sum(criteriaNumberOfBreakPoints)+2*numAlt)
  
  if (!is.null(alternativesRanks)){
    
    # determine now in which order the alternatives should be treated for the constraints
    indexOrder <- c()
    orderedAlternativesRanks <- sort(alternativesRanks)
    tmpRanks1 <- alternativesRanks
    tmpRanks2 <- alternativesRanks
    
    while (length(orderedAlternativesRanks) != 0){
      # search for the alternatives of lowest rank
      tmpIndex <- which(alternativesRanks == orderedAlternativesRanks[1])
      for (j in 1:length(tmpIndex))
        indexOrder<-c(indexOrder,tmpIndex[j])
      # remove the rank which has been dealt with now
      orderedAlternativesRanks<-orderedAlternativesRanks[-which(orderedAlternativesRanks==orderedAlternativesRanks[1])]
    }
    
    for (i in 1:(length(alternativesRanks)-1)){
      if (alternativesRanks[indexOrder[i]] == alternativesRanks[indexOrder[i+1]]){
        # then the alternatives are indifferent and their overall values are equal
        indifferenceConstraints <- rbind(indifferenceConstraints, a[indexOrder[i],] - a[indexOrder[i+1],])
      }
      else{
        # then the first alternative i is ranked better than the second one i+1 and i has an overall value higher than i+1
        preferenceConstraints <- rbind(preferenceConstraints, a[indexOrder[i],] - a[indexOrder[i+1],])
      } 
    }
  }
  
  if (!is.null(alternativesPreferences)){
    
    for (i in 1:dim(alternativesPreferences)[1]){
      preferenceConstraints <- rbind(preferenceConstraints, a[which(rownames(performanceTable)==alternativesPreferences[i,1]),] - a[which(rownames(performanceTable)==alternativesPreferences[i,2]),])
    }   
  }
  
  if (!is.null(alternativesIndifferences)){
    
    for (i in 1:dim(alternativesIndifferences)[1]){
      indifferenceConstraints <- rbind(indifferenceConstraints, a[which(rownames(performanceTable)==alternativesIndifferences[i,1]),] - a[which(rownames(performanceTable)==alternativesIndifferences[i,2]),])
    }   
  }
  
  # add this to the constraints matrix mat
  
  mat<-rbind(preferenceConstraints,indifferenceConstraints)
  
  # right hand side of this part of mat
  
  rhs <- c()
  
  if (dim(preferenceConstraints)[1]!=0){
    for (i in (1:dim(preferenceConstraints)[1]))
      rhs<-c(rhs,epsilon)
  }
  
  if (dim(indifferenceConstraints)[1]!=0){
    for (i in (1:dim(indifferenceConstraints)[1]))
      rhs<-c(rhs,0)
  }
  # direction of the inequality for this part of mat
  
  dir <- c()
  
  if (dim(preferenceConstraints)[1]!=0){
    for (i in (1:dim(preferenceConstraints)[1]))
      dir<-c(dir,">=")
  }
  
  if (dim(indifferenceConstraints)[1]!=0){
    for (i in (1:dim(indifferenceConstraints)[1]))
      dir<-c(dir,"==")
  }
  
  
  # -------------------------------------------------------
  
  # now the monotonicity constraints on the value functions
  
  monotonicityConstraints<-matrix(nrow=0, ncol=sum(criteriaNumberOfBreakPoints)+2*numAlt)
  
  for (i in 1:length(criteriaNumberOfBreakPoints)){
    for (j in 1:(criteriaNumberOfBreakPoints[i]-1)){
      tmp<-rep(0,sum(criteriaNumberOfBreakPoints)+2*numAlt)
      if (i==1)
        pos <- j
      else
        pos<-sum(criteriaNumberOfBreakPoints[1:(i-1)])+j
      tmp[pos] <- -1
      tmp[pos+1] <- 1
      monotonicityConstraints <- rbind(monotonicityConstraints, tmp)
    }
  }
  
  
  
  # add this to the constraints matrix mat
  
  mat<-rbind(mat,monotonicityConstraints)
  
  # the direction of the inequality
  
  for (i in (1:dim(monotonicityConstraints)[1]))
    dir<-c(dir,">=")
  
  # the right hand side of this part of mat
  
  for (i in (1:dim(monotonicityConstraints)[1]))
    rhs<-c(rhs,0)
  
  # -------------------------------------------------------
  
  # normalization constraint for the upper values of the value functions (sum = 1)
  
  tmp<-rep(0,sum(criteriaNumberOfBreakPoints)+2*numAlt)
  
  for (i in 1:length(criteriaNumberOfBreakPoints)){
    if (i==1)
      pos <- criteriaNumberOfBreakPoints[i]
    else
      pos<-sum(criteriaNumberOfBreakPoints[1:(i-1)])+criteriaNumberOfBreakPoints[i]
    tmp[pos] <- 1
  }
  
  # add this to the constraints matrix mat
  
  mat<-rbind(mat,tmp)
  
  # the direction of the inequality
  
  dir<-c(dir,"==")
  
  # the right hand side of this part of mat
  
  rhs<-c(rhs,1)
  
  # -------------------------------------------------------
  
  # now the normalizaiton constraints for the lower values of the value functions (= 0)
  
  minValueFunctionsConstraints<-matrix(nrow=0, ncol=sum(criteriaNumberOfBreakPoints)+2*numAlt)
  
  for (i in 1:length(criteriaNumberOfBreakPoints)){
    tmp<-rep(0,sum(criteriaNumberOfBreakPoints)+2*numAlt)
    if (i==1)
      pos <- i
    else
      pos<-sum(criteriaNumberOfBreakPoints[1:(i-1)])+1
    tmp[pos] <- 1
    minValueFunctionsConstraints <- rbind(minValueFunctionsConstraints,tmp)
  }
  
  # add this to the constraints matrix mat
  
  mat<-rbind(mat,minValueFunctionsConstraints)
  
  # the direction of the inequality
  
  for (i in (1:dim(minValueFunctionsConstraints)[1]))
    dir<-c(dir,"==")
  
  # the right hand side of this part of mat
  
  for (i in (1:dim(minValueFunctionsConstraints)[1]))
    rhs<-c(rhs,0)
  
  # -------------------------------------------------------
  
  lpSolution <- Rglpk_solve_LP(obj, mat, dir, rhs)
  
  # -------------------------------------------------------
  
  # create a structure containing the value functions
  
  valueFunctions <- list()
  
  for (i in 1:length(criteriaNumberOfBreakPoints)){
    tmp <- c() 
    if (i==1)
      pos <- 0
    else
      pos<-sum(criteriaNumberOfBreakPoints[1:(i-1)])
    for (j in 1:criteriaNumberOfBreakPoints[i]){
      tmp <- c(tmp,lpSolution$solution[pos+j])
    }
    tmp<-rbind(criteriaBreakPoints[[i]],tmp)
    colnames(tmp)<- NULL
    rownames(tmp) <- c("x","y")
    valueFunctions <- c(valueFunctions,list(tmp))
  }
  
  names(valueFunctions) <- colnames(performanceTable)
  
  # it might happen on certain computers that these value functions 
  # do NOT respect the monotonicity constraints (especially because of too small differences and computer arithmetics)
  # therefore we check if they do, and if not, we "correct" them
  
  for (i in 1:numCrit){
    for (j in 1:(criteriaNumberOfBreakPoints[i]-1)){
      if (valueFunctions[[i]][2,j] > valueFunctions[[i]][2,j+1]){
        valueFunctions[[i]][2,j+1] <- valueFunctions[[i]][2,j] 
      }
    }
  }
  
  
  # -------------------------------------------------------
  
  overallValues <- as.vector(t(a[,1:sum(criteriaNumberOfBreakPoints)]%*%lpSolution$solution[1:sum(criteriaNumberOfBreakPoints)]))
  
  names(overallValues) <- rownames(performanceTable)
  
  # -------------------------------------------------------
  
  # the error values for each alternative (sigma)
  
  errorValuesPlus <- as.vector(lpSolution$solution[(sum(criteriaNumberOfBreakPoints)+1):(sum(criteriaNumberOfBreakPoints)+numAlt)])
  errorValuesMinus <- as.vector(lpSolution$solution[(sum(criteriaNumberOfBreakPoints) + numAlt + 1):(sum(criteriaNumberOfBreakPoints)+2*numAlt)])
  
  names(errorValuesPlus) <- rownames(performanceTable)
  names(errorValuesMinus) <- rownames(performanceTable)
  
  # -------------------------------------------------------
  
  # the ranks of the alternatives 
  
  outRanks <- rank(-overallValues, ties.method="min")
  
  # -------------------------------------------------------
  
  if ((numAlt >= 3) && !is.null(alternativesRanks))
    tau = cor(alternativesRanks,outRanks,method="kendall")
  else
    tau = NULL
  
  # prepare the output
  
  out <- list(optimum = round(lpSolution$optimum, digits=5), valueFunctions = valueFunctions, overallValues = round(overallValues, digits=5), ranks = outRanks, errors = list(sigmaPlus = round(errorValuesPlus, digits=5), sigmaMinus = round(errorValuesMinus, digits=5)), Kendall = tau)
  
  
  # -------------------------------------------------------
  
  # post-optimality analysis if the optimum is found and if kPostOptimality is not NULL, i.e. the solution space is not empty
  
  minWeights <- NULL
  maxWeights <- NULL
  averageValueFunctions <- NULL
  
  
  if (!is.null(kPostOptimality) && (lpSolution$optimum == 0)){
    
    # add F \leq F* + k(F*) to the constraints, where F* is the optimum and k(F*) is a positive threshold, which is a small proportion of F*
    
    mat <- rbind(mat,obj)
    dir <- c(dir,"<=")
    rhs <- c(rhs,kPostOptimality)
    
    minWeights <- c()
    maxWeights <- c()
    combinedSolutions <- c()
    
    for (i in 1:numCrit){
      
      # first maximize the best ui for each criterion, then minimize it
      # this gives the interval of variation for each weight
      # the objective function : the first elements correspond to the ui's, the last one to the sigmas
      
      obj<-rep(0,sum(criteriaNumberOfBreakPoints))
      obj<-c(obj,rep(0,2*numAlt))
      
      if (i==1)
        pos <- criteriaNumberOfBreakPoints[i]
      else
        pos<-sum(criteriaNumberOfBreakPoints[1:(i-1)])+criteriaNumberOfBreakPoints[i]
      
      obj[pos] <- 1
      
      lpSolutionMin <- Rglpk_solve_LP(obj, mat, dir, rhs)
      lpSolutionMax <- Rglpk_solve_LP(obj, mat, dir, rhs, max=TRUE)
      
      minWeights <- c(minWeights,lpSolutionMin$optimum)
      maxWeights <- c(maxWeights,lpSolutionMax$optimum)
      combinedSolutions <- rbind(combinedSolutions,lpSolutionMin$solution)
      combinedSolutions <- rbind(combinedSolutions,lpSolutionMax$solution)
    }
    
    names(minWeights) <- colnames(performanceTable)
    names(maxWeights) <- colnames(performanceTable)
    
    # calculate the average value function, for which each component is the average value obtained for each of the programs above
    averageSolution <- apply(combinedSolutions,2,mean)
    
    # create a structure containing the average value functions
    
    averageValueFunctions <- list()
    
    for (i in 1:length(criteriaNumberOfBreakPoints)){
      tmp <- c() 
      if (i==1)
        pos <- 0
      else
        pos<-sum(criteriaNumberOfBreakPoints[1:(i-1)])
      for (j in 1:criteriaNumberOfBreakPoints[i]){
        tmp <- c(tmp,averageSolution[pos+j])
      }
      tmp<-rbind(criteriaBreakPoints[[i]],tmp)
      colnames(tmp)<- NULL
      rownames(tmp) <- c("x","y")
      averageValueFunctions <- c(averageValueFunctions,list(tmp))
    }
    
    names(averageValueFunctions) <- colnames(performanceTable)
    
  }
  
  out <- c(out, list(minimumWeightsPO = minWeights, maximumWeightsPO = maxWeights, averageValueFunctionsPO = averageValueFunctions))
  
  return(out)
}
