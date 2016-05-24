## Experimental research in evolutionary computation
## author: thomas.bartz-beielstein@fh-koeln.de
## http://www.springer.com/3-540-32026-1
##
## Copyright (C) 2004-2011  T. Bartz-Beielstein, C. Lasarczyk
## This program is free software;
## you can redistribute it and/or modify it under the terms of the GNU 
## General Public License as published by the Free Software Foundation;
## either version 3 of the License,
## or (at your option) any later version.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of 
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
## You should have received a copy of the GNU General Public License along 
## with this program; if not, see <http://www.gnu.org/licenses/>.
##
#### Verteilt das Budget optimal auf die Punkte, we consider minimization problems here
### spotOcba.R
### Parameters:
### samp.mean = vector of mean values, length nd
### samp.var = vector of variances, length nd
### samp.count = vector of repeats performed already, length nd
### budget.add = additional number of repeats, distributed among the nd design points 
### iz = indifference zone
### verbose = verbosity
#############################################

##################################################################################
#' Optimal Computing Budget Allocation OCBA for SPOT
#'
#' Spreads the budget in an optimal way for the different design points,
#' considering a minimization problem
#'
#' @param samp.mean vector of mean values, length nd
#' @param samp.var vector of variances, length nd
#' @param samp.count vector of repeats performed already, length nd
#' @param budget.add additional number of repeats, distributed among the nd design points 
#' @param iz indifference zone
#' @param verbose verbosity {0|1|2|3} 0 is no printing
#' 
#' @seealso  \code{\link{spotGenerateSequentialDesign}} \code{\link{spotRepeatsOcba}}
#' @export
###################################################################################
spotOcba <- function(samp.mean, samp.var, samp.count, budget.add, iz=NA, verbose=0){
	spotWriteLines(verbose,2,"  Entering spotOcba")
	nd <- length(samp.mean)
	### create a double-precision vector of the specified length nd with each element equal to 0:
	ratio <- numeric(nd);
	## Nach dieser Runde verbrauchtes Budget
	budget.total <- sum(samp.count) + budget.add
	samp.order <- order(samp.mean)
	## Es kann mehrere beste Settings geben, dann wird das naechst schlechteste
	## zum Zweitbesten gemacht, damit die aus dem Unterschied berechnete
	## Indiffenzzone nicht Null wird. 
	best.mean <- min(samp.mean)
	nrOfBest <- sum(samp.mean==best.mean)
	if (nrOfBest > 1) warning("More than a single best setting")
	### Now we consider the case that every design point has the same function value:
	if (nrOfBest == nd){
		samp.additional <- numeric(nd)
		samp.additional[1] <- budget.add
		return(samp.additional)
	}	
	best.id <- samp.order[1:nrOfBest]
	second.id <- samp.order[nrOfBest+1]
	ratio[second.id] <- 1.0
	## Fuer all ausser die ersten beiden
	if (verbose>1) 
		print(best.id)
	### we consider not the best or second best:
	for(i in samp.order[(nrOfBest + 2):nd]) {
		if (is.na(iz)) {
			temp <- (samp.mean[second.id] - best.mean) /
				(samp.mean[i] - best.mean) 
		} else {
			temp <- samp.mean[second.id] - best.mean
			temp <- max(temp, iz)
			temp1 <- samp.mean[i] - best.mean
			temp1 <- max(temp1, iz)
			temp <- temp/temp1
		}
		ratio[i] <- temp^2*samp.var[i]/samp.var[second.id]
	}
	if (verbose>1) print(ratio)
	if (is.na(iz)) {
		temp <- sum(ratio * ratio / samp.var)
		ratio[best.id] <- sqrt(samp.var[best.id] * temp)
	} else {
		ratio[best.id] <- samp.var[best.id] / samp.var[second.id]
	}
	morerun <- rep(TRUE,nd)
	moreAlloc <- TRUE
	budget.temp <- budget.total
	samp.desired <- numeric(nd)
	while(moreAlloc) {
		moreAlloc <- FALSE
		# print(ratio)
		ratio.sum <- sum(ratio[morerun])
		if (verbose>1)
			print(ratio.sum)
		samp.desired[morerun] <- floor((budget.temp/ratio.sum)*ratio[morerun])
		# print(morerun);
		## Der Sachverhalt "Mehr als Genug" wird nur einmal festgestellt, da dann ein Anpassung stattfindet
		samp.moreThanEnough <- (samp.desired < samp.count)
		samp.desired[samp.moreThanEnough] <- samp.count[samp.moreThanEnough]
		moreAlloc <- any(samp.moreThanEnough)
		morerun[samp.moreThanEnough] <- FALSE
		if (moreAlloc) budget.temp <- budget.total - sum(samp.desired[!morerun])
	}
	samp.additional <- samp.desired - samp.count
	## Differenz wird dem Besten zugeschlagen
	budget.temp <- sum(samp.desired)
	samp.additional[best.id] <- samp.additional[best.id] + (budget.total - budget.temp)
	spotWriteLines(verbose,2,"  Leaving spotOcba")
	samp.additional
}