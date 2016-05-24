#############################################################################
#
# Copyright Patrick Meyer, Sébastien Bigaret and Richard Hodgett, 2015
#
# Contributors:
#   Patrick Meyer <patrick.meyer@telecom-bretagne.eu>
#   Sebastien Bigaret <sebastien.bigaret@telecom-bretagne.eu>
#   Richard Hodgett <r.e.hodgett@leeds.ac.uk>
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

TOPSIS <- function(performanceTable, criteriaWeights, criteriaMinMax, positiveIdealSolutions = NULL, negativeIdealSolutions = NULL, alternativesIDs = NULL, criteriaIDs = NULL){
	
	## check the input data

        if (!(is.null(alternativesIDs) || is.vector(alternativesIDs)))
        	stop("alternatives IDs should be in a vector")
        	
        if (!(is.null(criteriaIDs) || is.vector(criteriaIDs)))
        	stop("criteria IDs should be in a vector")

	if (!((is.matrix(performanceTable) || (is.data.frame(performanceTable))))) 
        stop("performanceTable must be a matrix or a data frame")

	if (!(length(criteriaWeights) == ncol(performanceTable))) 
		stop("the number of criteriaWeights must equal the number of columns in the performanceTable")

	if (missing(criteriaMinMax))
		stop("the input criteriaMinMax is required.")

	## filter the performance table and the criteria according to the given alternatives and criteria
	
	if (!is.null(alternativesIDs)) performanceTable <- performanceTable[alternativesIDs,]
	
	if (!is.null(criteriaIDs)) 	{
							performanceTable <- performanceTable[,criteriaIDs]
							criteriaWeights <- criteriaWeights[criteriaIDs]
							if (!missing(positiveIdealSolutions)) positiveIdealSolutions <- positiveIdealSolutions[criteriaIDs]
							if (!missing(negativeIdealSolutions)) negativeIdealSolutions <- negativeIdealSolutions[criteriaIDs]
							criteriaMinMax <- criteriaMinMax[criteriaIDs]
						}
	
	critno <- length(criteriaWeights)
	altno <- nrow(performanceTable)

	## Calculate the weighted normalised matrix
	divby <- c(1:critno)
	for (i in 1:critno)
		{
			divby[i] <- sqrt(sum(performanceTable[,i]^2))
		}
	normalisedm <- t(t(performanceTable) / divby)
	wnm <- t(t(normalisedm) * criteriaWeights)

	## Identify positive and negative ideal solutions
	pis <- c(1:critno)
	nis <- c(1:critno)
	if (missing(positiveIdealSolutions) || missing(negativeIdealSolutions))
		{
		for (i in 1:critno)
			{
				if (criteriaMinMax[i] == "max")
				{
					pis[i] <- max(wnm[,i])
					nis[i] <- min(wnm[,i])
				}
				else
				{
					pis[i] <- min(wnm[,i])
					nis[i] <- max(wnm[,i])
				}
			}
		}
	else
		{
		## check the input data is correct
		if (!(length(positiveIdealSolutions) == length(negativeIdealSolutions) || length(positiveIdealSolutions) == critno)) 
		stop("the number of postive and negaitve ideal solutions need to equal the number of alternaitves.")
		pis <- positiveIdealSolutions
		nis <- negativeIdealSolutions
		}

	## Identify separation from positive and negative ideal solutions
	spis <- sweep(wnm,2,pis)^2
	snis <- sweep(wnm,2,nis)^2	
	spisv <- c(1:altno)
	snisv <- c(1:altno)

	for (i in 1:altno)
			{
				spisv[i] <- sqrt(sum(spis[i,]))
				snisv[i] <- sqrt(sum(snis[i,]))
			}

	## Calculate results
	results <- c(1:altno)
	for (i in 1:altno)
			{
				results[i] <- snisv[i] / (snisv[i] + spisv[i])
			}
	names(results) <- rownames(performanceTable)
		return(results)
}