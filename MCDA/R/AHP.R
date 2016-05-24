#############################################################################
#
# Copyright Patrick Meyer, Sabastien Bigaret and Richard Hodgett, 2015
#
# Contributors:
#   Patrick Meyer <patrick.meyer@telecom-bretagne.eu>
#   Sabastien Bigaret <sebastien.bigaret@telecom-bretagne.eu>
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

AHP <- function(criteriaWeightsPairwiseComparisons, alternativesPairwiseComparisonsList){
	
	## check the input data is correct

	if (!((is.matrix(criteriaWeightsPairwiseComparisons) || (is.data.frame(criteriaWeightsPairwiseComparisons))))) 
        stop("criteriaWeightsPairwiseComparisons must be a matrix or a data frame")

	if (!(nrow(criteriaWeightsPairwiseComparisons) == ncol(criteriaWeightsPairwiseComparisons))) 
		stop("criteriaWeightsPairwiseComparisons must be a square matrix or a data frame")

	if(!all(criteriaWeightsPairwiseComparisons == t(1/criteriaWeightsPairwiseComparisons)))
		stop("criteriaWeightsPairwiseComparisons must be a reciprocal matrix (i.e. value on one side must = 1/value)")

	if (!(length(alternativesPairwiseComparisonsList) >= 2)) 
        stop("list alternativesPairwiseComparisonsList must contain at least 2 matrices or data frames")

	size <- 0;
	for (i in 1:length(alternativesPairwiseComparisonsList))
		{
			if (!((is.matrix(alternativesPairwiseComparisonsList[i][[1]]) || (is.data.frame(alternativesPairwiseComparisonsList[i][[1]]))))) 
        		stop("all elements in the list alternativesPairwiseComparisonsList must be a matrix or a data frame")

			if (!(nrow(alternativesPairwiseComparisonsList[i][[1]]) == ncol(alternativesPairwiseComparisonsList[i][[1]]))) 
			stop("all elements in the list alternativesPairwiseComparisonsList must be a square matrix or a data frame")

			if(!all(alternativesPairwiseComparisonsList[i][[1]] == t(1/alternativesPairwiseComparisonsList[i][[1]])))
			stop("all elements in the list scoresmatrixlist must be a reciprocal matrix (i.e. value on one side must = 1/value)")
			
			if (i == 1)
			{
				size <- nrow(alternativesPairwiseComparisonsList[i][[1]])
			}
			else
			{
				if (!(nrow(alternativesPairwiseComparisonsList[i][[1]]) == size)) 
				stop("all matrices or data frames in list alternativesPairwiseComparisonsList must be the same size.")
			}		
		}

	critno <- nrow(criteriaWeightsPairwiseComparisons)
	altno <- nrow(alternativesPairwiseComparisonsList[1][[1]])

	## Estimate the principle eigenvector of the weights matrix to 10 digits
	
	pairwisematrix <-  criteriaWeightsPairwiseComparisons  %*% criteriaWeightsPairwiseComparisons
	sumrows <- rowSums(criteriaWeightsPairwiseComparisons)
	sumtotal <- sum(sumrows)
	normalisedsumrows <- sumrows / sumtotal	
	previous <- vector()
	
	while (!identical(round(previous, digits = 10),round(normalisedsumrows, digits = 10)))
					{	
	previous <- normalisedsumrows
	pairwisematrix <-  pairwisematrix  %*% pairwisematrix
	sumrows <- rowSums(pairwisematrix)
	sumtotal <- sum(sumrows)
	normalisedsumrows <- sumrows / sumtotal
					}
	weights <- normalisedsumrows

	## Estimate the principle eigenvectors of each of the score matrices to 10 digits

	savedscores <- matrix(nrow=critno,ncol=altno)
	for (i in 1:length(alternativesPairwiseComparisonsList))
		{
			pairwisematrix <-  alternativesPairwiseComparisonsList[i][[1]]  %*% alternativesPairwiseComparisonsList[i][[1]]
			sumrows <- rowSums(alternativesPairwiseComparisonsList[i][[1]])
			sumtotal <- sum(sumrows)
			normalisedsumrows <- sumrows / sumtotal	
			previous <- vector()
	
			while (!identical(round(previous, digits = 10),round(normalisedsumrows, digits = 10)))
					{	
				previous <- normalisedsumrows
				pairwisematrix <-  pairwisematrix  %*% pairwisematrix
				sumrows <- rowSums(pairwisematrix)
				sumtotal <- sum(sumrows)
				normalisedsumrows <- sumrows / sumtotal
					}
			
			savedscores[i,] <- normalisedsumrows
		}

	## Calculate the results

	results <- matrix(nrow=1,ncol=altno)
	
	for (i in 1:altno)
		{
		altscore <- 0
	  	for (j in 1:critno)
			{
				altscore <- altscore + (weights[j] * savedscores[j,i])
			}
		results[1,i] <- altscore
		}
  results<-as.vector(results)
  names(results)<-row.names(alternativesPairwiseComparisonsList[[1]])
	return(results)
}