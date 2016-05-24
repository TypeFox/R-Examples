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

plotMRSortSortingProblem <- function(performanceTable, categoriesLowerProfiles, assignments, criteriaMinMax, criteriaUBs, criteriaLBs, alternativesIDs = NULL, criteriaIDs = NULL){
  
  ## check the input data
  
  if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
    stop("wrong performanceTable, should be a matrix or a data frame")
  
  if (!(is.matrix(categoriesLowerProfiles)))
    stop("categoriesLowerProfiles should be a matrix")
  
  if (!(is.vector(assignments)))
    stop("assignments should be a vector")
  
  if (!(is.vector(criteriaMinMax)))
    stop("criteriaMinMax should be a vector")
  
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
    categoriesLowerProfiles <- categoriesLowerProfiles[,criteriaIDs]
    criteriaUBs <- criteriaUBs[criteriaIDs]
    criteriaLBs <- criteriaLBs[criteriaIDs]
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
  
  normalizedPerformanceTable <- matrix(nrow=numAlt,ncol=numCrit)
  
  for (j in 1:numAlt){
    for (i in 1:numCrit){
      if(criteriaMinMax[i] == "min")
        normalizedPerformanceTable[j,i] <- 1-(performanceTable[j,i]-criteriaLBs[i])/(criteriaUBs[i]-criteriaLBs[i])
      else
        normalizedPerformanceTable[j,i] <- (performanceTable[j,i]-criteriaLBs[i])/(criteriaUBs[i]-criteriaLBs[i])
    }
  }
  
  normalizedProfiles <- matrix(nrow=numCat,ncol=numCrit)
  
  for (j in 1:numCat){
    for (i in 1:numCrit){
      if(criteriaMinMax[i] == "min")
        normalizedProfiles[j,i] <- 1-(categoriesLowerProfiles[j,i]-criteriaLBs[i])/(criteriaUBs[i]-criteriaLBs[i])
      else
        normalizedProfiles[j,i] <- (categoriesLowerProfiles[j,i]-criteriaLBs[i])/(criteriaUBs[i]-criteriaLBs[i])
    }
  }
  
  col.rainbow <- rainbow(numAlt)
  
  palette(col.rainbow)
  
  ylim=c(min(normalizedPerformanceTable)-0.1,max(normalizedPerformanceTable)+0.1)
  
  # Expand right side of clipping rect to make room for the legend
  # par(xpd=T, mar=par()$mar+c(0,0,0,6))
  
  # par(mfrow=c(1,2))
  
  #layout(c(1, 2), widths=c(7, 1))
  
  layout(matrix(c(1,2),1), widths = c(4,2))
  
  plot(1:numCrit,normalizedProfiles[1,],
       type="l",
       col="black",
       ylim=ylim,
       xlab = "",
       ylab="",
       xaxt="n",
       yaxt="n",
       lwd=2,
       lty=2)
  # the other profiles (except the lower one)
  for (i in (2:numCat-1))
    lines(1:numCrit,normalizedProfiles[i,],col="black", lwd=2, lty=2)
  for (i in 1:numCrit){
    lines(c(i,i),ylim, col="gray")
  }
  axis(1,at=c(1:numCrit),labels=colnames(performanceTable))
  
  for (i in 1:(numCat-1)){
    text(c(1:numCrit), normalizedProfiles[i,], labels = categoriesLowerProfiles[i,], pos=1)
  }
  
  text(1, ylim[2], labels = dimnames(performanceTable)[[2]][1], pos=4)
  
   for (i in 2:(numCrit)){
    text(i, ylim[2], labels = dimnames(performanceTable)[[2]][i], pos=2)
   }
    
  for (i in (1:numAlt))
    points(1:numCrit,normalizedPerformanceTable[i,],type="b",pch=which(assignments[rownames(performanceTable)[i]]==rownames(categoriesLowerProfiles)), col=i)
  
  par(mar=c(0, 0, 0, 0))
  plot.new()
  
  legend("center", c(rownames(performanceTable),rownames(categoriesLowerProfiles)), cex=0.8, col=c(c(1:numAlt),rep("black",numCat)), 
         lwd=1, bty="n",pch=c(rep(0,numAlt),c(1:numCat)));
  
}
