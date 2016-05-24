# 		QStat: tool for calculating the Q statistic of Goeman's global test for metabolomic pathways, part of the 'MetStaT' package  
#		Copyright (C) 2012 Diana Hendrickx and Tim Dorscheidt
#		
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation, either version 3 of the License, or
#		(at your option) any later version.
#		
#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.
#		
#		You should have received a copy of the GNU General Public License
#		along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# 		Email: g.zwanenburg@uva.nl ('MetStaT' contact person) or tdorscheidt@gmail.com
###############################################################################

# Wrapper function that calculates the Q statistic and performs a permutation test.
QStat.Calculate <- function(X, y.boolean, permutations = "all") {
	
	# Internal core Qstat calculation method. Adapted with permission from the original Matlab code of the author. See referenced article in copyright notice for further details.
	QStat.GetQValue <- function(X, y.boolean) {
		q <- sum(y.boolean)
		t <- length(y.boolean)
		mu <- q / t
		m <- dim(X)[2]
		Z <- y.boolean - mu
		R <- X%*%t(X)/m
		Q <- t(Z)%*%R%*%Z/(mu*(1-mu))
		Q
	}
	
	result <- list()
	if (class(y.boolean)!="logical") {
		y.boolean <- MetStaT.ConvertToNumericClasses(y.boolean,new.classes=c(1,0))==1
	}
	result$y.boolean <- y.boolean
	result$X <- X
	
	Q <- QStat.GetQValue(X, y.boolean)
	result$Q <- Q
	
	if (permutations!=0) {
		result$p <- QStat.PermutationTest(result, no.permutations = permutations)
	} else {
		result$p <- NA
	}
	
	result
}

# Method for performing a permutation test on an already performed Q statistic.
QStat.PermutationTest <- function(Qresult, no.permutations = "all", quietly = FALSE) {
	# Internal method for calculating the number of combinations possible in a binary vector of n length with a subset of k non-default values.  
	QStat.CombinationsPossible <- function(n, k) {
		gamma(n+1)/(gamma(k+1)*gamma(n-k+1)) # gamma is the most accessible R factorial(x) method, but it excludes x itself in the product series, therefore x needs to be increased by 1
	}
	# Internal core Qstat calculation method. Adapted with permission from the original Matlab code of the author. See referenced article in copyright notice for further details.
	QStat.GetQValue <- function(X, y.boolean) {
		q <- sum(y.boolean)
		t <- length(y.boolean)
		mu <- q / t
		m <- dim(X)[2]
		Z <- y.boolean - mu
		R <- X%*%t(X)/m
		Q <- t(Z)%*%R%*%Z/(mu*(1-mu))
		Q
	}
	# confirm which method of going through permutations needs to be used
	max.permutations <- QStat.CombinationsPossible(length(Qresult$y.boolean),sum(Qresult$y.boolean))
	if (no.permutations!="all" && is.numeric(no.permutations) && no.permutations > max.permutations) {
		if (no.permutations==0) return(NA)
		# OPTION 1: randomly select a fixed number of permutations to be used for the test
		count.higher.Qs <- 0 # keep track of how many of the permutated Q calculations are above the non-permutated one
		if (!quietly) print(paste("Performing exhaustive test of all",max.permutations,"permutations."))
		for (p in 1:no.permutations) {
			y.perm <- sample(Qresult$y.boolean)
			Q.perm <- QStat.GetQValue(Qresult$X, y.perm)
			if (Q.perm >= Qresult$Q) {
				count.higher.Qs <- count.higher.Qs + 1		
			}
		}
		return(count.higher.Qs/no.permutations)
	} else { # OPTION 2: exhaustively go through all permutation options for the test
		if (!quietly) print(paste("Performing test using",no.permutations,"random permutations."))
		# define a recursive function that iterates through all possible permutations (couple of orders of magnitude faster than the combn function by R)
		QStat.NextBooleanCombination <- function(location.of.trues, max, which.true.to.relocate = 1) {
			if (which.true.to.relocate > length(location.of.trues)) {return(NULL)} # no further combinations possible
			if (location.of.trues[which.true.to.relocate] <= max - which.true.to.relocate) { # can I move the current 'true' in question?
				location.of.trues[which.true.to.relocate] <- location.of.trues[which.true.to.relocate] + 1 # if so, then move it
				if (which.true.to.relocate > 1) { # move all 'trues' that are further ahead as close to the current one as possible
					for (c in (which.true.to.relocate-1):1) {
						location.of.trues[c] <- location.of.trues[c+1]+1 
					}
				}
			} else {
				return(QStat.NextBooleanCombination(location.of.trues, max, which.true.to.relocate + 1)) # try moving a 'true' earlier in the vector
			}
			location.of.trues
		}
		y.perm.empty <- rep(FALSE,length(Qresult$y.boolean)) # initialize an empty boolean y vector 
		location.of.trues <- c(sum(Qresult$y.boolean):1) # initialize an array that contains the default starting positions of the trues in the y vector
		count.higher.Qs <- 0 # keep track of how many of the permutated Q calculations are above the non-permutated one 
		# now loop through all possible permutations
		while(!is.null(location.of.trues)) {
			y.perm <- y.perm.empty
			y.perm[location.of.trues]=TRUE
			Q.perm <- QStat.GetQValue(Qresult$X, y.perm)
			if (Q.perm >= Qresult$Q) {
				count.higher.Qs <- count.higher.Qs + 1		
			}
			location.of.trues <- QStat.NextBooleanCombination(location.of.trues, length(Qresult$y.boolean))
		}
		return(count.higher.Qs/max.permutations)
	}
}
