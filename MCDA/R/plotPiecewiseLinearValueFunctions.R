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

plotPiecewiseLinearValueFunctions <- function(valueFunctions, criteriaIDs = NULL){
	
	## check the input data
  
	if (!(is.list(valueFunctions)))
        	stop("valueFunctions should be a list")
  
	if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
	  stop("criteriaIDs should be a vector")
  
	## filter the data according to the given criteria
	
	if (!is.null(criteriaIDs)){
	  valueFunctions <- valueFunctions[criteriaIDs]
	}
  
	if (is.null(valueFunctions[[1]]))
	  stop("no value functions left to plot")
  else 
    numCrit <- length(valueFunctions)
  
	# plotting symbol and color 
  
	par(pch=22, col="red")
  
	# determine how many plots per row and column
  
  if (numCrit <= 2)
    par(mfrow=c(1,2)) 
	else
	  par(mfrow=c(ceiling(log2(numCrit)),ceiling(log2(numCrit)))) 
  
  # plot the functions
    
	for(i in 1:numCrit){
	  heading = names(valueFunctions)[i]
	  plot(valueFunctions[[i]]["x",], valueFunctions[[i]]["y",], type="n", main=heading,xlab="", ylab="") 
	  lines(valueFunctions[[i]]["x",], valueFunctions[[i]]["y",], type="b") 
	}
  
}
