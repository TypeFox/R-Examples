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

UTADIS <- function(performanceTable, criteriaMinMax, criteriaNumberOfBreakPoints, alternativesAssignments, categoriesRanks, epsilon, criteriaLBs=NULL, criteriaUBs=NULL, alternativesIDs = NULL, criteriaIDs = NULL, categoriesIDs = NULL){
  
  ## check the input data
  
  if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
    stop("wrong performanceTable, should be a matrix or a data frame")
  
  if (!(is.null(alternativesAssignments) || is.vector(alternativesAssignments)))
    stop("alternativesRanks should be a vector")
  
  if (!(is.null(categoriesRanks) || is.vector(categoriesRanks)))
    stop("categoriesRanks should be a vector")
  
  if (!(is.vector(criteriaMinMax)))
    stop("criteriaMinMax should be a vector")
  
  if (!(is.vector(criteriaNumberOfBreakPoints)))
    stop("criteriaNumberOfBreakPoints should be a vector")
  
  if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
    stop("alternativesIDs should be in a vector")
  
  if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
    stop("criteriaIDs should be in a vector")
  
  if (!(is.null(categoriesIDs) || is.vector(categoriesIDs)))
    stop("categoriesIDs should be in a vector")
  
  if (!(is.null(criteriaLBs) || is.vector(criteriaLBs)))
    stop("criteriaLBs should be in a vector")
  
  if (!(is.null(criteriaUBs) || is.vector(criteriaUBs)))
    stop("criteriaUBs should be in a vector")
  
  ## filter the data according to the given alternatives and criteria
  ## in alternativesIDs and criteriaIDs
  
  if (!is.null(alternativesIDs)){
    performanceTable <- performanceTable[alternativesIDs,] 
    alternativesAssignments <- alternativesAssignments[alternativesIDs]
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
  
  ## filter the data according to the given cateogries
  ## in categoriesIDs
  ## also update the performanceTable to remove alternatives
  ## which have no assignments anymore
  
  if (!is.null(categoriesIDs)){
    tmp<-lapply(alternativesAssignments, function(x) x == categoriesIDs)
    tmp2<-c()
    for (i in 1:length(tmp)){
      tmp2<-c(tmp2,any(tmp[[i]]))
    }
    alternativesAssignments<-alternativesAssignments[tmp2]
    categoriesRanks <- categoriesRanks[categoriesRanks=categoriesIDs]
    performanceTable <- performanceTable[names(alternativesAssignments),] 
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
  
  # if there are no assignments examples left
  # we stop here
  
  if (length(alternativesAssignments)==0)
    stop("after filtering alternativesAssignments is empty")
  
  # if there are no categories ranks left
  # we stop here
  
  if (length(categoriesRanks) == 0)
    stop("after filtering categoriesRanks is empty")
  
  # check if categoriesRanks and alternativesAssignments are compatible
  # i.e. if for each assignment example, the category has a rank
  
  if (length(setdiff(unique(alternativesAssignments),names(categoriesRanks))) != 0)
    stop("alternativesAssignments contains categories which have no rank in categoriesRanks")
  
  
  # -------------------------------------------------------
  
  numCrit <- dim(performanceTable)[2]
  
  numAlt <- dim(performanceTable)[1]
  
  # define how many categories we have and how they are named (after all the filtering processes)
  
  categoriesIDs <- unique(alternativesAssignments)
  
  numCat <- length(categoriesIDs)
  
  if (numCat <= 1)
    stop("only 1 or less category left after filtering")
  
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
  
  # the objective function : the first elements correspond to the ui's, then the sigmas, and to finish, we have the category thresholds (one lower threshold per category, none for the lowest category)
  
  obj<-rep(0,sum(criteriaNumberOfBreakPoints))
  
  obj<-c(obj,rep(1,2*numAlt))
  
  obj<-c(obj,rep(0,numCat-1))
  
  # -------------------------------------------------------
  
  assignmentConstraintsLBs <- matrix(nrow=0, ncol=sum(criteriaNumberOfBreakPoints)+2*numAlt + numCat - 1)
  assignmentConstraintsUBs <- matrix(nrow=0, ncol=sum(criteriaNumberOfBreakPoints)+2*numAlt + numCat - 1)
  
  # for each assignment example, write its constraint
  
  # categoriesRanks might contain ranks which are higher than 1,2,3, ... (for example if some categories have been filtered out)
  # that's why we recalculate the ranks of the filtered categories here
  # that will give us the positions in the lower bounds vector of the categories
  
  newCategoriesRanks <- rank(categoriesRanks)
  
  for (i in 1:length(alternativesAssignments)){
    # determine which lower bound should be activated for the comparison
    # none if it is an assignment in the lowest category
    if (newCategoriesRanks[alternativesAssignments[i]] != max(newCategoriesRanks)){
      tmp <- rep(0,numCat-1)
      tmp[newCategoriesRanks[alternativesAssignments[i]]] <- -1
      assignmentConstraintsLBs<-rbind(assignmentConstraintsLBs, c(a[which(rownames(performanceTable) == names(alternativesAssignments)[i]),],tmp))
    }
    if (newCategoriesRanks[alternativesAssignments[i]] != min(newCategoriesRanks)){
      tmp <- rep(0,numCat-1)
      tmp[newCategoriesRanks[alternativesAssignments[i]]-1] <- -1
      assignmentConstraintsUBs<-rbind(assignmentConstraintsUBs, c(a[which(rownames(performanceTable) == names(alternativesAssignments)[i]),],tmp))
    }
  }
  
  
  # add this to the constraints matrix mat
  
  mat<-rbind(assignmentConstraintsLBs,assignmentConstraintsUBs)
  
  # right hand side of this part of mat
  
  rhs <- c()
  
  if (dim(assignmentConstraintsLBs)[1]!=0){
    for (i in (1:dim(assignmentConstraintsLBs)[1]))
      rhs<-c(rhs,0)
  }
  
  if (dim(assignmentConstraintsUBs)[1]!=0){
    for (i in (1:dim(assignmentConstraintsUBs)[1]))
      rhs<-c(rhs,-epsilon)
  }
  # direction of the inequality for this part of mat
  
  dir <- c()
  
  if (dim(assignmentConstraintsLBs)[1]!=0){
    for (i in (1:dim(assignmentConstraintsLBs)[1]))
      dir<-c(dir,">=")
  }
  
  if (dim(assignmentConstraintsUBs)[1]!=0){
    for (i in (1:dim(assignmentConstraintsUBs)[1]))
      dir<-c(dir,"<=")
  }
  
  
  # -------------------------------------------------------
  
  # now the monotonicity constraints on the value functions
  
  monotonicityConstraints<-matrix(nrow=0, ncol=sum(criteriaNumberOfBreakPoints)+2*numAlt + numCat - 1)
  
  for (i in 1:length(criteriaNumberOfBreakPoints)){
    for (j in 1:(criteriaNumberOfBreakPoints[i]-1)){
      tmp<-rep(0,sum(criteriaNumberOfBreakPoints)+2*numAlt + numCat - 1)
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
  
  tmp<-rep(0,sum(criteriaNumberOfBreakPoints)+2*numAlt + numCat - 1)
  
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
  
  minValueFunctionsConstraints<-matrix(nrow=0, ncol=sum(criteriaNumberOfBreakPoints)+2*numAlt + numCat - 1)
  
  for (i in 1:length(criteriaNumberOfBreakPoints)){
    tmp<-rep(0,sum(criteriaNumberOfBreakPoints)+2*numAlt + numCat - 1)
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
  
  # the categories' lower bounds
  
  categoriesLBs <- as.vector(lpSolution$solution[(sum(criteriaNumberOfBreakPoints)+2*numAlt+1):(sum(criteriaNumberOfBreakPoints)+2*numAlt+numCat-1)])
  names(categoriesLBs) <- names(newCategoriesRanks[1:numCat-1])
  
  #   
  #   # -------------------------------------------------------
  #   
  #   # the ranks of the alternatives 
  #   
  #   outRanks <- rank(-overallValues, ties.method="min")
  #   
  #   # -------------------------------------------------------
  #   
  #   if ((numAlt >= 3) && !is.null(alternativesRanks))
  #     tau = Kendall(alternativesRanks,outRanks)$tau[1]
  #   else
  #     tau = NULL
  #   
  #   # prepare the output
  #   
  out <- list(optimum = round(lpSolution$optimum, digits=5), valueFunctions = valueFunctions, overallValues = round(overallValues, digits=5), categoriesLBs = categoriesLBs, errors = list(sigmaPlus = round(errorValuesPlus, digits=5), sigmaMinus = round(errorValuesMinus, digits=5)))
  #   
  #   
  #   # -------------------------------------------------------
  #   
  #   # post-optimality analysis if the optimum is found and if kPostOptimality is not NULL, i.e. the solution space is not empty
  #   
  #   minWeights <- NULL
  #   maxWeights <- NULL
  #   averageValueFunctions <- NULL
  #   
  #   
  #   if (!is.null(kPostOptimality) && (lpSolution$optimum == 0)){
  #     
  #     # add F \leq F* + k(F*) to the constraints, where F* is the optimum and k(F*) is a positive threshold, which is a small proportion of F*
  #     
  #     mat <- rbind(mat,obj)
  #     dir <- c(dir,"<=")
  #     rhs <- c(rhs,kPostOptimality)
  #     
  #     minWeights <- c()
  #     maxWeights <- c()
  #     combinedSolutions <- c()
  #     
  #     for (i in 1:numCrit){
  #       
  #       # first maximize the best ui for each criterion, then minimize it
  #       # this gives the interval of variation for each weight
  #       # the objective function : the first elements correspond to the ui's, the last one to the sigmas
  #       
  #       obj<-rep(0,sum(criteriaNumberOfBreakPoints))
  #       obj<-c(obj,rep(0,2*numAlt))
  #       
  #       if (i==1)
  #         pos <- criteriaNumberOfBreakPoints[i]
  #       else
  #         pos<-sum(criteriaNumberOfBreakPoints[1:(i-1)])+criteriaNumberOfBreakPoints[i]
  #       
  #       obj[pos] <- 1
  #       
  #       lpSolutionMin <- Rglpk_solve_LP(obj, mat, dir, rhs)
  #       lpSolutionMax <- Rglpk_solve_LP(obj, mat, dir, rhs, max=TRUE)
  #       
  #       minWeights <- c(minWeights,lpSolutionMin$optimum)
  #       maxWeights <- c(maxWeights,lpSolutionMax$optimum)
  #       combinedSolutions <- rbind(combinedSolutions,lpSolutionMin$solution)
  #       combinedSolutions <- rbind(combinedSolutions,lpSolutionMax$solution)
  #     }
  #     
  #     names(minWeights) <- colnames(performanceTable)
  #     names(maxWeights) <- colnames(performanceTable)
  #     
  #     # calculate the average value function, for which each component is the average value obtained for each of the programs above
  #     averageSolution <- apply(combinedSolutions,2,mean)
  #     
  #     # create a structure containing the average value functions
  #     
  #     averageValueFunctions <- list()
  #     
  #     for (i in 1:length(criteriaNumberOfBreakPoints)){
  #       tmp <- c() 
  #       if (i==1)
  #         pos <- 0
  #       else
  #         pos<-sum(criteriaNumberOfBreakPoints[1:(i-1)])
  #       for (j in 1:criteriaNumberOfBreakPoints[i]){
  #         tmp <- c(tmp,averageSolution[pos+j])
  #       }
  #       tmp<-rbind(criteriaBreakPoints[[i]],tmp)
  #       colnames(tmp)<- NULL
  #       rownames(tmp) <- c("x","y")
  #       averageValueFunctions <- c(averageValueFunctions,list(tmp))
  #     }
  #     
  #     names(averageValueFunctions) <- colnames(performanceTable)
  #     
  #   }
  #   
  #   out <- c(out, list(minimumWeightsPO = minWeights, maximumWeightsPO = maxWeights, averageValueFunctionsPO = averageValueFunctions))
  
  return(out)
}
