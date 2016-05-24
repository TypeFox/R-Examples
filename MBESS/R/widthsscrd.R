# Find the variance of the conditions difference with or without a covariate
.varCRDDiff <- function(nclus, nindiv, prtreat, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL) {
	# Fill in tauy or sigmay by the specification of totalvar and iccy
	if(is.null(tauy)) {
		if(is.null(iccy)) {
			tauy <- totalvar - sigma2y
		} else {
			tauy <- iccy * totalvar
		}
	}
	if(is.null(sigma2y)) {
		if(is.null(iccy)) {
			sigma2y <- totalvar - tauy
		} else {
			sigma2y <- (1 - iccy) * totalvar
		}	
	}
	if(length(nindiv) > 1) {
		sigma2y <- sigma2y * (1 - r2within)
		tauy <- tauy * (1 - r2between)
		D <- 1/sigma2y
		F <- 0 - (tauy)/(sigma2y * (sigma2y + nindiv*tauy))
		ntreat <- round(prtreat*nclus)
		prtreat <- ntreat/nclus
		terms <- nindiv * (D + nindiv*F)
		termstotal <- sum(terms)
		termstreatment <- sum(terms[1:ntreat])
		vargamma1 <- termstotal/(termstotal*termstreatment - termstreatment^2)
		varclusmean <- vargamma1 * (nclus * prtreat * (1 - prtreat))
	} else {
		ntreat <- round(prtreat*nclus)
		prtreat <- ntreat/nclus
		varclusmean <- ((sigma2y * (1 - r2within)) / nindiv) + (tauy * (1 - r2between)) # variance of a single predicted cluster mean
	}
	if(!is.null(assurance)) { # Adjust the variance of a single predicted cluster mean when the degree of assurance is specified
		df <- nclus - 2 - numpredictor
		varclusmean <- varclusmean * qchisq(assurance, df) / df
	}
	denominator <- nclus * prtreat * (1 - prtreat) # The other terms in the variance of the conditions difference formula
	return(varclusmean/denominator)
} 

# Find the total cost of the whole cluster randomized design
.costCRD <- function(nclus, nindiv, cluscost, indivcost, diffsize = NULL) {
	nindiv[nindiv == "> 100000"] <- Inf # If the number of individuals is '> 100000', the price is set as infinity
	nindiv <- as.numeric(nindiv)
	if(!is.null(diffsize)) nindiv <- .findNindivVec(nindiv, diffsize, nclus)
	result <- NULL
	if(length(nindiv) > 1) {
		result <- (nclus * cluscost) + sum(nindiv * indivcost)
	} else {
		result <- nclus * (cluscost + (nindiv * indivcost))
	}
	return(result)
}

# Find width of the CI of the unstandardized conditions difference
.findWidthCRDDiff <- function(nclus, nindiv, prtreat, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize = NULL) {
	nclus <- as.numeric(nclus)
	prtreat <- round(nclus * prtreat)/nclus
	nindiv <- as.numeric(nindiv)
	if(!is.null(diffsize)) nindiv <- .findNindivVec(nindiv, diffsize, nclus)
	v <- .varCRDDiff(nclus=nclus, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance)
	alpha <- 1 - conf.level
	width <- 2 * sqrt(v) * qt(1 - alpha/2, nclus - 2 - numpredictor)
	return(width)
}

.findNindivVec <- function(nindiv, diffsize, nclus) {
	isInteger <- all(round(diffsize) == diffsize)
	result <- NULL
	if(isInteger) {
		result <- rep(nindiv + diffsize, length.out=nclus)
	} else {
		result <- rep(round(nindiv * diffsize), length.out=nclus)
	}
	result[result < 1] <- 1
	result
}

# Find the number of clusters given the specified width and the cluster size
.findNclusCRDDiff <- function(width, nindiv, prtreat, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize=NULL) {
	nclus <- seq(100, 1000, 100) # Starting values
	result <- sapply(nclus, .findWidthCRDDiff, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
	
	# Find the position where the width posited
	if(all(width > result)) { # If the specified width is larger than all calculated widths, the number of clusters should be less than 100
		nclus <- seq(ceiling(2 * (1/min(prtreat, 1 - prtreat))) + numpredictor, 100, 1) # The minimum is to make sure that there are at least two groups in each condition
		result <- sapply(nclus, .findWidthCRDDiff, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
		if(all(width > result)) {
			return(3 + numpredictor)
		} else if (all(width < result)) {
			return(100)
		} else {
			return(nclus[which(width > result)[1]])
		}	
	} else if (all(width < result)) { # If the specified width is smaller than all calculated width, the number of clusters should be larger than 1000
		start <- 1000
		repeat {
			nclus <- seq(start, start + 1000, 100)
			result <- sapply(nclus, .findWidthCRDDiff, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
			if(all(width > result)) {
				return(start)	
			} else if (all(width < result)) {
				start <- start + 1000
			} else {
				minval <- nclus[which(width > result)[1] - 1]
				maxval <- nclus[which(width > result)[1]]
				nclus <- seq(minval, maxval, 1)
				result <- sapply(nclus, .findWidthCRDDiff, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
				if(all(width > result)) {
					return(minval)
				} else if (all(width < result)) {
					return(maxval)
				} else {
					return(nclus[which(width > result)[1]])
				}		
			}
		}
	} else { # The specified width is somewhere in between the calculated width
		minval <- nclus[which(width > result)[1] - 1]
		maxval <- nclus[which(width > result)[1]]
		nclus <- seq(minval, maxval, 1)
		result <- sapply(nclus, .findWidthCRDDiff, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
		if(all(width > result)) {
			return(minval)
		} else if (all(width < result)) {
			return(maxval)
		} else {
			return(nclus[which(width > result)[1]])
		}			
	}
}

# Find the cluster size given the specified width and the number of clusters
.findNindivCRDDiff <- function(width, nclus, prtreat, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize = NULL) {
	baremaximum <- .findWidthCRDDiff(nclus=nclus, nindiv=100000, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
	if(width < baremaximum) return("> 100000")
	nindiv <- seq(100, 1000, 100)
	result <- sapply(nindiv, .findWidthCRDDiff, nclus=nclus, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
	if(all(width > result)) {
		nindiv <- seq(2, 100, 1)
		result <- sapply(nindiv, .findWidthCRDDiff, nclus=nclus, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
		if(all(width > result)) {
			return(2)
		} else if (all(width < result)) {
			return(100)
		} else {
			return(nindiv[which(width > result)[1]])
		}	
	} else if (all(width < result)) {
		start <- 1000
		repeat {
			nindiv <- seq(start, start + 1000, 100)
			result <- sapply(nindiv, .findWidthCRDDiff, nclus=nclus, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
			if(all(width > result)) {
				return(start)	
			} else if (all(width < result)) {
				start <- start + 1000
			} else {
				minval <- nindiv[which(width > result)[1] - 1]
				maxval <- nindiv[which(width > result)[1]]
				nindiv <- seq(minval, maxval, 1)
				result <- sapply(nindiv, .findWidthCRDDiff, nclus=nclus, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
				if(all(width > result)) {
					return(minval)
				} else if (all(width < result)) {
					return(maxval)
				} else {
					return(nindiv[which(width > result)[1]])
				}		
			}
		}
	} else {
		minval <- nindiv[which(width > result)[1] - 1]
		maxval <- nindiv[which(width > result)[1]]
		nindiv <- seq(minval, maxval, 1)
		result <- sapply(nindiv, .findWidthCRDDiff, nclus=nclus, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
		if(all(width > result)) {
			return(minval)
		} else if (all(width < result)) {
			return(maxval)
		} else {
			return(nindiv[which(width > result)[1]])
		}			
	}
}

# Find the least expensive combination of the number of clusters and cluster size given the specified width
.findMinCostCRDDiff <- function(width, cluscost=0, indivcost=1, prtreat, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize = NULL) {
	nclus <- seq(100, 1100, 100)
	repeat {
		nindiv <- sapply(nclus, .findNindivCRDDiff, width=width, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
		cost <- mapply(.costCRD, nclus=nclus, nindiv=nindiv, MoreArgs=list(cluscost=cluscost, indivcost=indivcost, diffsize=diffsize), SIMPLIFY=TRUE)
		if(!all(cost == Inf)) break
		nclus <- nclus + 1000
	}
	posmin <- which(cost == min(cost))
	if(length(posmin) > 1) posmin <- posmin[1]
	if(posmin == 1) {
		nclus <- seq(ceiling(2 * (1/min(prtreat, 1 - prtreat))) + numpredictor, 200, 1)
		nindiv <- sapply(nclus, .findNindivCRDDiff, width=width, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
		cost <- mapply(.costCRD, nclus=nclus, nindiv=nindiv, MoreArgs=list(cluscost=cluscost, indivcost=indivcost, diffsize=diffsize), SIMPLIFY=TRUE)
		posmin <- which(cost == min(cost))
		if(length(posmin) > 1) posmin <- posmin[1]
		return(c(nclus[posmin], nindiv[posmin], cost[posmin]))	
	} else if (posmin == length(cost)) {
		start <- nclus[length(nclus) - 1]
		repeat {
			nclus <- seq(start, start + 1100, 100)
			nindiv <- sapply(nclus, .findNindivCRDDiff, width=width, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
			cost <- mapply(.costCRD, nclus=nclus, nindiv=nindiv, MoreArgs=list(cluscost=cluscost, indivcost=indivcost, diffsize=diffsize), SIMPLIFY=TRUE)
			posmin <- which(cost == min(cost))
			if(length(posmin) > 1) posmin <- posmin[1]
			if(posmin == 1) {
				return(c(nclus[posmin], nindiv[posmin], cost[posmin]))	
			} else if (posmin == length(cost)) {
				start <- start + 1000
			} else {
				minval <- nclus[posmin - 1]
				maxval <- nclus[posmin + 1]
				nclus <- seq(minval, maxval, 1)
				nindiv <- sapply(nclus, .findNindivCRDDiff, width=width, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
				cost <- mapply(.costCRD, nclus=nclus, nindiv=nindiv, MoreArgs=list(cluscost=cluscost, indivcost=indivcost, diffsize=diffsize), SIMPLIFY=TRUE)
				posmin <- which(cost == min(cost))
				if(length(posmin) > 1) posmin <- posmin[1]
				return(c(nclus[posmin], nindiv[posmin], cost[posmin]))			
			}
		}
	} else {
		minval <- nclus[posmin - 1]
		maxval <- nclus[posmin + 1]
		nclus <- seq(minval, maxval, 1)
		nindiv <- sapply(nclus, .findNindivCRDDiff, width=width, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
		cost <- mapply(.costCRD, nclus=nclus, nindiv=nindiv, MoreArgs=list(cluscost=cluscost, indivcost=indivcost, diffsize=diffsize), SIMPLIFY=TRUE)
		posmin <- which(cost == min(cost))
		if(length(posmin) > 1) posmin <- posmin[1]
		return(c(nclus[posmin], nindiv[posmin], cost[posmin]))			
	}
}

# Find the cluster size given the budget and the number of clusters
.findNindivCRDBudget <- function(budget, nclus, cluscost, indivcost, diffsize = NULL) {
	nindiv <- seq(100, 1000, 100)
	result <- sapply(nindiv, .costCRD, nclus=nclus, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	if(all(budget < result)) {
		nindiv <- seq(2, 100, 1)
		result <- sapply(nindiv, .costCRD, nclus=nclus, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
		if(all(budget < result)) {
			return(NA)
		} else if (all(budget > result)) {
			return(100)
		} else {
			index <- which(budget >= result)
			return(nindiv[index[length(index)]])
		}	
	} else if (all(budget > result)) {
		start <- 1000
		repeat {
			nindiv <- seq(start, start + 1000, 100)
			result <- sapply(nindiv, .costCRD, nclus=nclus, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
			if(all(budget < result)) {
				return(start)	
			} else if (all(budget > result)) {
				start <- start + 1000
			} else {
				minval <- nindiv[which(budget < result)[1] - 1]
				maxval <- nindiv[which(budget < result)[1]]
				nindiv <- seq(minval, maxval, 1)
				result <- sapply(nindiv, .costCRD, nclus=nclus, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
				if(all(budget < result)) {
					return(minval)
				} else if (all(budget > result)) {
					return(maxval)
				} else {
					index <- which(budget >= result)
					return(nindiv[index[length(index)]])
				}		
			}
		}
	} else {
		minval <- nindiv[which(budget < result)[1] - 1]
		maxval <- nindiv[which(budget < result)[1]]
		nindiv <- seq(minval, maxval, 1)
		result <- sapply(nindiv, .costCRD, nclus=nclus, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
		if(all(budget < result)) {
			return(minval)
		} else if (all(budget > result)) {
			return(maxval)
		} else {
			index <- which(budget >= result)
			return(nindiv[index[length(index)]])
		}			
	}
}

# Find the least width combination of the number of clusters and cluster size given the budget
.findMinWidthCRDDiff <- function(budget, cluscost=0, indivcost=1, prtreat, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize = NULL) {
	FUN <- function(nclus, nindiv) {
		.findWidthCRDDiff(nclus=nclus, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level=conf.level, diffsize=diffsize)
	}
	nclus <- seq(100, 1100, 100)
	nindiv <- sapply(nclus, .findNindivCRDBudget, budget=budget, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	if(all(is.na(nindiv))) {
		posmin <- 1
	} else {
		resultwidth <- rep(NA, length(nclus))
		resultwidth[!is.na(nindiv)] <- mapply(FUN, nclus=nclus[!is.na(nindiv)], nindiv=nindiv[!is.na(nindiv)], SIMPLIFY=TRUE)
		posmin <- which(resultwidth == min(resultwidth, na.rm=TRUE))
	}
	if(posmin == 1) {
		nclus <- seq(ceiling(2 * (1/min(prtreat, 1 - prtreat))) + numpredictor, 200, 1)
		nindiv <- sapply(nclus, .findNindivCRDBudget, budget=budget, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
		resultwidth <- rep(NA, length(nclus))
		resultwidth[!is.na(nindiv)] <- mapply(FUN, nclus=nclus[!is.na(nindiv)], nindiv=nindiv[!is.na(nindiv)])
		posmin <- which(resultwidth == min(resultwidth, na.rm=TRUE))
		return(c(nclus[posmin], nindiv[posmin], resultwidth[posmin]))	
	} else if (posmin == length(resultwidth)) {
		start <- nclus[length(nclus) - 1]
		repeat {
			nclus <- seq(start, start + 1100, 100)
			nindiv <- sapply(nclus, .findNindivCRDBudget, budget=budget, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
			resultwidth <- rep(NA, length(nclus))
			resultwidth[!is.na(nindiv)] <- mapply(FUN, nclus=nclus[!is.na(nindiv)], nindiv=nindiv[!is.na(nindiv)])
			posmin <- which(resultwidth == min(resultwidth, na.rm=TRUE))
			if(posmin == 1) {
				return(c(nclus[posmin], nindiv[posmin], resultwidth[posmin]))	
			} else if (posmin == length(resultwidth)) {
				start <- start + 1000
			} else {
				minval <- nclus[posmin - 1]
				maxval <- nclus[posmin + 1]
				nclus <- seq(minval, maxval, 1)
				nindiv <- sapply(nclus, .findNindivCRDBudget, budget=budget, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
				resultwidth <- rep(NA, length(nclus))
				resultwidth[!is.na(nindiv)] <- mapply(FUN, nclus=nclus[!is.na(nindiv)], nindiv=nindiv[!is.na(nindiv)])
				posmin <- which(resultwidth == min(resultwidth, na.rm=TRUE))
				return(c(nclus[posmin], nindiv[posmin], resultwidth[posmin]))			
			}
		}
	} else {
		minval <- nclus[posmin - 1]
		maxval <- nclus[posmin + 1]
		nclus <- seq(minval, maxval, 1)
		nindiv <- sapply(nclus, .findNindivCRDBudget, budget=budget, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
		resultwidth <- rep(NA, length(nclus))
		resultwidth[!is.na(nindiv)] <- mapply(FUN, nclus=nclus[!is.na(nindiv)], nindiv=nindiv[!is.na(nindiv)])
		posmin <- which(resultwidth == min(resultwidth, na.rm=TRUE))
		return(c(nclus[posmin], nindiv[posmin], resultwidth[posmin]))			
	}
}

# Find the number of clusters given the budget and the cluster size
.findNclusCRDBudget <- function(budget, nindiv, cluscost, indivcost, diffsize = NULL) {
	nclus <- seq(100, 1000, 100)
	result <- sapply(nclus, .costCRD, nindiv=nindiv, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	if(all(budget < result)) {
		nclus <- seq(4, 100, 1)
		result <- sapply(nclus, .costCRD, nindiv=nindiv, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
		if(all(budget < result)) {
			return(NA)
		} else if (all(budget > result)) {
			return(100)
		} else {
			index <- which(budget >= result)
			return(nclus[index[length(index)]])
		}	
	} else if (all(budget > result)) {
		start <- 1000
		repeat {
			nclus <- seq(start, start + 1000, 100)
			result <- sapply(nclus, .costCRD, nindiv=nindiv, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
			if(all(budget < result)) {
				return(start)	
			} else if (all(budget > result)) {
				start <- start + 1000
			} else {
				minval <- nclus[which(budget < result)[1] - 1]
				maxval <- nclus[which(budget < result)[1]]
				nclus <- seq(minval, maxval, 1)
				result <- sapply(nclus, .costCRD, nindiv=nindiv, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
				if(all(budget < result)) {
					return(minval)
				} else if (all(budget > result)) {
					return(maxval)
				} else {
					index <- which(budget >= result)
					return(nclus[index[length(index)]])
				}		
			}
		}
	} else {
		minval <- nclus[which(budget < result)[1] - 1]
		maxval <- nclus[which(budget < result)[1]]
		nclus <- seq(minval, maxval, 1)
		result <- sapply(nclus, .costCRD, nindiv=nindiv, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
		if(all(budget < result)) {
			return(minval)
		} else if (all(budget > result)) {
			return(maxval)
		} else {
			index <- which(budget >= result)
			return(nclus[index[length(index)]])
		}			
	}
}


# Get a nice report for a public function
.reportCRD <- function(nclus, nindiv, width=NULL, cost=NULL, es=FALSE, estype=0, assurance=NULL, diffsize = NULL) {
	cat(paste("The number of clusters is ", nclus, ".\n", sep=""))
	if(is.null(diffsize)) {
		cat(paste("The cluster size is ", nindiv, ".\n", sep=""))
	} else {
		cat("The cluster size is as follows:\n")
		out <- data.frame(table(.findNindivVec(as.numeric(nindiv), diffsize, nclus)))
		colnames(out) <- c("Cluster Size", "Frequency")
		print(out)
	}
	if(!is.null(width)) {
		if(es) {
			eslab <- NULL
			if(estype == 0) {
				eslab <- "total"
			} else if (estype == 1) {
				eslab <- "individual-level"
			} else {
				eslab <- "cluster-level"
			}
			if(is.null(assurance)) {
				cat(paste("The expected width of ", eslab, " effect size is ", round(width, 4), ".\n", sep=""))
			} else {
				cat(paste("The width of ", eslab, " effect size with ", round(assurance, 2), " assurance is ", round(width, 4), ".\n", sep=""))
			}
		} else {
			if(is.null(assurance)) {
				cat(paste("The expected width of unstandardized conditions difference is ", round(width, 4), ".\n", sep=""))
			} else {
				cat(paste("The width of unstandardized conditions difference with ", round(assurance, 2), " assurance is ", round(width, 4), ".\n", sep=""))	
			}
		}
	}
	if(!is.null(cost)) cat(paste("The budget is ", cost, ".\n", sep=""))
}


ss.aipe.crd.nclus.fixedwidth <- function(width, nindiv, prtreat, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, cluscost=NULL, indivcost=NULL, diffsize = NULL) {
	nclus <- .findNclusCRDDiff(width=width, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	calculatedWidth <- .findWidthCRDDiff(nclus=nclus, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	calculatedCost <- NULL
	if(!is.null(cluscost) && !is.null(indivcost)) calculatedCost <- .costCRD(nclus, nindiv, cluscost, indivcost, diffsize = diffsize)
	.reportCRD(nclus, nindiv, calculatedWidth, cost=calculatedCost, es=FALSE, estype=0, assurance=assurance, diffsize = diffsize) 
	invisible(nclus)
}


ss.aipe.crd.nindiv.fixedwidth <- function(width, nclus, prtreat, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, cluscost=NULL, indivcost=NULL, diffsize = NULL) {
	nindiv <- .findNindivCRDDiff(width=width, nclus=nclus, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	if(nindiv == "> 100000") stop("With the current number of clusters, it is impossible to achieve the target width. Please increase the number of clusters.")
	calculatedWidth <- .findWidthCRDDiff(nclus=nclus, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	calculatedCost <- NULL
	if(!is.null(cluscost) && !is.null(indivcost)) calculatedCost <- .costCRD(nclus, nindiv, cluscost, indivcost, diffsize = diffsize)
	.reportCRD(nclus, nindiv, calculatedWidth, cost=calculatedCost, es=FALSE, estype=0, assurance=assurance, diffsize = diffsize) 
	invisible(nindiv)
}


ss.aipe.crd.nclus.fixedbudget <- function(budget, nindiv, cluscost = 0, indivcost = 1, prtreat = NULL, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize = NULL) {
	nclus <- .findNclusCRDBudget(budget=budget, nindiv=nindiv, cluscost=cluscost, indivcost=indivcost, diffsize = diffsize)
	calculatedWidth <- NULL
	if(!is.null(prtreat)) {
		calculatedWidth <- .findWidthCRDDiff(nclus=nclus, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	}
	calculatedCost <- .costCRD(nclus, nindiv, cluscost=cluscost, indivcost=indivcost, diffsize = diffsize)
	.reportCRD(nclus, nindiv, calculatedWidth, cost=calculatedCost, es=FALSE, estype=0, assurance=assurance, diffsize = diffsize) 
	invisible(nclus)
}


ss.aipe.crd.nindiv.fixedbudget <- function(budget, nclus, cluscost = 0, indivcost = 1, prtreat = NULL, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize = NULL) {
	nindiv <- .findNindivCRDBudget(budget=budget, nclus=nclus, cluscost=cluscost, indivcost=indivcost, diffsize = diffsize)
	calculatedWidth <- NULL
	if(!is.null(prtreat)) {
		calculatedWidth <- .findWidthCRDDiff(nclus=nclus, nindiv=nindiv, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	}
	calculatedCost <- .costCRD(nclus, nindiv, cluscost=cluscost, indivcost=indivcost, diffsize = diffsize)
	.reportCRD(nclus, nindiv, calculatedWidth, cost=calculatedCost, es=FALSE, estype=0, assurance=assurance, diffsize = diffsize) 
	invisible(nindiv)
}


ss.aipe.crd.both.fixedbudget <- function(budget, cluscost=0, indivcost=1, prtreat, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize = NULL) {
	result <- .findMinWidthCRDDiff(budget=budget, cluscost=cluscost, indivcost=indivcost, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	calculatedCost <- .costCRD(result[1], result[2], cluscost=cluscost, indivcost=indivcost, diffsize = diffsize)
	.reportCRD(result[1], result[2], result[3], cost=calculatedCost, es=FALSE, estype=0, assurance=assurance, diffsize = diffsize) 
	invisible(result[1:2])
}


ss.aipe.crd.both.fixedwidth <- function(width, cluscost=0, indivcost=1, prtreat, tauy=NULL, sigma2y=NULL, totalvar=NULL, iccy=NULL, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, diffsize = NULL) {
	result <- .findMinCostCRDDiff(width=width, cluscost=cluscost, indivcost=indivcost, prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	calculatedWidth <- .findWidthCRDDiff(nclus=result[1], nindiv=result[2], prtreat=prtreat, tauy=tauy, sigma2y=sigma2y, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	.reportCRD(result[1], result[2], calculatedWidth, cost=result[3], es=FALSE, estype=0, assurance=assurance, diffsize = diffsize) 
	invisible(result[1:2])
}

# Create dataset from the CRD model
.createDataCRD <- function(nclus, ntreatclus, nindiv, iccy, es, estype = 1, totalvar=1, covariate=FALSE, iccz=NULL, r2within=NULL, r2between=NULL, totalvarz = 1, diffsize = NULL) {
if(!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")

	if(!is.null(diffsize)) { 
		nindiv <- .findNindivVec(nindiv, diffsize, nclus)
	} else {
		nindiv <- rep(nindiv, each=nclus)
	}
	id <- do.call(c, mapply(rep, 1:nclus, nindiv, SIMPLIFY=FALSE)) #rep(1:nclus, each=nindiv) # id variable
	x <- c(rep(1, ntreatclus), rep(0, nclus - ntreatclus)) # Create a grouping variable
	tau <- iccy * totalvar # Between-group variance
	sigma <- (1 - iccy) * totalvar # Within-group variance
	gamma1 <- NULL # Unstandardized variance
	if(estype == 0) { # Create unstandardized difference between conditions
		gamma1 <- es * sqrt(totalvar)
	} else if (estype == 1) {
		gamma1 <- es * sqrt(sigma)
	} else if (estype == 2) {
		gamma1 <- es * sqrt(tau)
	} else {
		stop("'estype' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
	}
	if(covariate) {
		if(iccz == 0 && r2between != 0) stop("Because the covariate varies at the level 1 only, the r-square at level 2 must be 0.")
		if(iccz == 1 && r2within != 0) stop("Because the covariate varies at the level 2 only, the r-square at level 1 must be 0.")
		tauz <- totalvarz * iccz # between-group variance of the covariate
		sigmaz <- totalvarz * (1 - iccz) # within-group variance of the covariate
		gammazb <- sqrt(r2between * tau / tauz) # the effect of the covariate on the between level
		gammazw <- sqrt(r2within * sigma / sigmaz) # the effect of the covariate on the within level
		tau <- (1 - r2between) * tau # Update residual between-group variance
		sigma <- (1 - r2within) * sigma # Update residual within-group variance
	}
	ybetween <- (gamma1 * x) + rnorm(nclus, 0, sqrt(tau)) # Create the between-level dependent variable
	if(covariate && iccz != 0) {
		zbetween <- rnorm(nclus, 0, sqrt(tauz)) # Create the between-level covariate
		ybetween <- ybetween + gammazb * zbetween # Add the covariate into the dependent variable
	}
	ybetween <- do.call(c, mapply(rep, ybetween, nindiv, SIMPLIFY=FALSE)) #rep(ybetween, each=nindiv) # Repeat the between-level dependent variable with nindiv times
	y <- ybetween + rnorm(sum(nindiv), 0, sqrt(sigma)) # Add the within-level dv
	if(covariate && iccz != 1) {
		zwithin <- rnorm(sum(nindiv), 0, sqrt(sigmaz)) # Create the within-level covariate
		y <- y + gammazw * zwithin # Add the covariate into the dv
	}
	x <- do.call(c, mapply(rep, x, nindiv, SIMPLIFY=FALSE)) #rep(x, each=nindiv) # Treatment group variables
	dat <- data.frame(id, y, x) # Attach the data frame
	if(covariate) {
		z <- NULL # Create a covariate
		if(iccz == 0) {
			z <- zwithin
		} else if(iccz == 1) {
			z <- do.call(c, mapply(rep, zbetween, nindiv, SIMPLIFY=FALSE)) #rep(zbetween, each = nindiv)
		} else {
			z <- do.call(c, mapply(rep, zbetween, nindiv, SIMPLIFY=FALSE)) + zwithin
		}
		if(iccz == 0) {
			zb <- 0
			zw <- z # If the intraclass correlation of z = 0, the within-level covariate is treated as a covariate
		} else {
			#zlist <- as.list(as.data.frame(matrix(z, nrow=nindiv))) # Make the group-mean centering on the covariate
			zlist <- split(z, id)
			zb <- do.call(c, mapply(rep, sapply(zlist, mean), nindiv, SIMPLIFY=FALSE)) #rep(sapply(zlist, mean), each=nindiv) 
			zw <- do.call(c, lapply(zlist, scale, scale=FALSE))
		}
		z <- data.frame(zw=zw, zb=zb)
		dat <- data.frame(dat, z) # Attach the covariate
	}
	rownames(dat) <- NULL
	return(dat)
}

# Change the data to wide-format
.wideFormat <- function(data, betweencol, withincol, idcol) {
	temp <- split(data[,withincol], data[,idcol]) # Separate the variables in the level 1 into a list which each element means level-2 units
	temp <- lapply(temp, function(x) as.vector(as.matrix(x))) # Transform matrix into a vector
	dataw <- do.call(rbind, temp) # Combine them together
	datab <- as.matrix(data[match(unique(data[,idcol]), data[,idcol]), betweencol]) # Concatenate the transformed data to the between-level variables
	colnames(datab) <- colnames(data)[betweencol] 
	nindiv <- nrow(data) / nrow(datab)
	colnames(dataw) <- paste(rep(colnames(data)[withincol], each=nindiv), rep(1:nindiv, length(withincol)), sep="") 
	return(data.frame(datab, dataw))
}

# Change the data to wide-format
.wideFormatUnequal <- function(data, betweencol, withincol, idcol) {
	temp <- split(data[,withincol], data[,idcol]) # Separate the variables in the level 1 into a list which each element means level-2 units
	# maxsize <- max(table(data[,idcol]))
	# FUN <- function(dat, size) {
		# row <- nrow(dat)
		# col <- ncol(dat)
		# dat <- as.matrix(dat)
		# return(rbind(dat, matrix(NA, size - row, col)))
	# }
	#temp <- lapply(temp, FUN, size=maxsize)
	temp <- lapply(temp, function(x) as.vector(as.matrix(x))) # Transform matrix into a vector
	
	size <- sapply(temp, length)/length(withincol)
	dataw <- lapply(split(temp, size), function(x) do.call(rbind, x)) # Combine them together
	datab <- split(data[match(unique(data[,idcol]), data[,idcol]), betweencol], size)
	resultdat <- mapply(data.frame, datab, dataw)
	varnamesw <- lapply(sapply(dataw, ncol)/length(withincol), function(x) paste(rep(colnames(data)[withincol], each=x), rep(1:x, length(withincol)), sep=""))
	varnames <- lapply(varnamesw, function(x) c(colnames(data)[betweencol], x))
	resultdat <- mapply(function(x, y) { colnames(x) <- y; x}, x = resultdat, y = varnames)
	return(resultdat)
}

# Create dataset from the CRD model and transform it to the wide format
.createDataCRDWide <- function(nclus, ntreatclus, nindiv, iccy, es, estype = 1, totalvar=1, covariate=FALSE, iccz=NULL, r2within=NULL, r2between=NULL, totalvarz = 1, diffsize = NULL) {
	dat <- .createDataCRD(nclus=nclus, ntreatclus=ntreatclus, nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=covariate, iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = totalvarz, diffsize = diffsize) # Long-format data
	datawide <- NULL # Create wide format data
	if(!is.null(diffsize)) {
		if(covariate) {
			if(iccz == 0) {
				datawide <- .wideFormatUnequal(dat, 3, c(2, 4), 1)
			} else if(iccz == 1) {
				datawide <- .wideFormatUnequal(dat, c(3, 5), 2, 1)
			} else {
				datawide <- .wideFormatUnequal(dat, c(3, 5), c(2, 4), 1)
			}
		} else {
			datawide <- .wideFormatUnequal(dat, 3, 2, 1)
		}
	} else {
		if(covariate) {
			if(iccz == 0) {
				datawide <- .wideFormat(dat, 3, c(2, 4), 1)
			} else if(iccz == 1) {
				datawide <- .wideFormat(dat, c(3, 5), 2, 1)
			} else {
				datawide <- .wideFormat(dat, c(3, 5), c(2, 4), 1)
			}
		} else {
			datawide <- .wideFormat(dat, 3, 2, 1)
		}	
	}
	return(datawide)
}

# Calculate a likelihood-based CI of ES
.likCIESCRD <- function(datawide, ylab, xlab, zwlab=NULL, zblab=NULL, estype=1, iccy=0.25, es=0.5, totalvar=1, covariate=FALSE, iccz=0.25, r2within=0.5, r2between=0.5, totalvarz = 1, conf.level = 0.95) {
	tau <- iccy * totalvar # Between-group variance
	sigma <- (1 - iccy) * totalvar # Within-group variance
	gamma1 <- NULL # Unstandardized conditions difference
	if(estype == 0) { # Create the unstandardized conditions difference based on the type of effect size
		gamma1 <- es * sqrt(totalvar)
	} else if (estype == 1) {
		gamma1 <- es * sqrt(sigma)
	} else if (estype == 2) {
		gamma1 <- es * sqrt(tau)
	} else {
		stop("'estype' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
	}
	if(covariate) {
		if(iccz == 0 && r2between != 0) stop("Because the covariate varies at the level 1 only, the r-square at level 2 must be 0.")
		if(iccz == 1 && r2within != 0) stop("Because the covariate varies at the level 2 only, the r-square at level 1 must be 0.")
		tauz <- totalvarz * iccz
		sigmaz <- totalvarz * (1 - iccz)
		gammazb <- sqrt(r2between * tau / tauz)
		gammazw <- sqrt(r2within * sigma / sigmaz)
		tau <- (1 - r2between) * tau # Update residual between-group variance
		sigma <- (1 - r2within) * sigma # Update residual within-group variance
	}
	probx <- sum(datawide[,xlab])/nrow(datawide) # Probability of the treatment condition

	if(!requireNamespace("OpenMx", quietly = TRUE)) stop("The package 'OpenMx' is needed; please install the package and try again.")
	
	
	latentlab <- c("intcept", "slope") # Intercept = Between-group variation of dv; slope = The within-level effect of the covariate
	if(is.null(zwlab)) latentlab <- "intcept"

	# Get the dimensions of matrix in the RAM notation
	frowlab <- c(ylab, xlab, zblab) # Manifest variables labels
	fcollab <- c(frowlab, latentlab) # Manifest+Latent variables labels
	lenrow <- length(frowlab) 
	lencol <- length(fcollab)
	
	# Regression coefficients matrix
	Alab <- matrix(NA, lencol, lencol)
	Aval <- matrix(0, lencol, lencol)
	Afree <- matrix(FALSE, lencol, lencol)
	colnames(Alab) <- colnames(Aval) <- colnames(Afree) <- rownames(Alab) <- rownames(Aval) <- rownames(Afree) <- fcollab
	
	Alab["intcept", xlab] <- "groupdiff" # unstandardized conditions difference
	if(!is.null(zblab)) Alab["intcept", zblab] <- "zbeffect" # The between-level effect of the covariate
	if(!is.null(zwlab)) Alab[ylab, "slope"] <- paste("data.", zwlab, sep="") # Definition variables putting the within-level covariate values in the factor loadings from the "slope" factor
	
	Aval["intcept", xlab] <- gamma1 
	if(!is.null(zblab)) Aval["intcept", zblab] <- gammazb 
	Aval[ylab, latentlab] <- 1
	
	Afree["intcept", xlab] <- TRUE
	if(!is.null(zblab)) Afree["intcept", zblab] <- TRUE

	# Covariance matrix
	Slab <- matrix(NA, lencol, lencol)
	Sval <- matrix(0, lencol, lencol)
	Sfree <- matrix(FALSE, lencol, lencol)
	colnames(Slab) <- colnames(Sval) <- colnames(Sfree) <- rownames(Slab) <- rownames(Sval) <- rownames(Sfree) <- fcollab
	diag(Slab)[1:length(ylab)] <- "l1error" # Resiual variance in the within level
	diag(Sval)[1:length(ylab)] <- sigma 
	diag(Sfree)[1:length(ylab)] <- TRUE
	Slab["intcept", "intcept"] <- "l2error" # Residual variance in the between level
	Sval["intcept", "intcept"] <- tau 
	Sfree["intcept", "intcept"] <- TRUE
	Slab[xlab, xlab] <- "varx" # The variance of the grouping variables
	Sval[xlab, xlab] <- probx * (1 - probx) 
	Sfree[xlab, xlab] <- TRUE
	if(!is.null(zblab)) {
		Slab[c(xlab, zblab), c(xlab, zblab)] <- "covxzb" # The covariance between between-level covariate and grouping variable
		Slab[xlab, xlab] <- "varx" # The variance of the grouping variable
		Slab[zblab, zblab] <- "varzb" # The variance of the between-level covariate
		Sval[c(xlab, zblab), c(xlab, zblab)] <- 0
		Sval[xlab, xlab] <- probx * (1 - probx) 
		Sval[zblab, zblab] <- tauz 
		Sfree[c(xlab, zblab), c(xlab, zblab)] <- TRUE
	}

	# The selection matrix
	Fval <- cbind(diag(lenrow), matrix(0, lenrow, length(latentlab)))
	Flab <- matrix(NA, lenrow, lencol)
	Ffree <- matrix(FALSE, lenrow, lencol)
	colnames(Flab) <- colnames(Fval) <- colnames(Ffree) <- fcollab
	rownames(Flab) <- rownames(Fval) <- rownames(Ffree) <- frowlab

	# The intercept vectors
	Mlab <- c(rep(NA, length(ylab)), "meanX") # The proportion of treatment group
	Mval <- c(rep(0, length(ylab)), probx) 
	Mfree <- c(rep(FALSE, length(ylab)), TRUE)
	
	if(!is.null(zblab)) {
		Mlab <- c(Mlab, "meanzb") # The average of the between-level covariate
		Mval <- c(Mval, 0)
		Mfree <- c(Mfree, TRUE)	
	}
	Mlab <- c(Mlab, "meanctrl") # The average of the control condition
	Mval <- c(Mval, 0)
	Mfree <- c(Mfree, TRUE)
	if(!is.null(zwlab)) {
		Mlab <- c(Mlab, "zweffect") # The average of the within-level covariate
		Mval <- c(Mval, gammazw) 
		Mfree <- c(Mfree, TRUE)	
	} 	
	
	Mlab <- matrix(Mlab, 1, lencol) # Change a vector to an one-row matrix
	Mval <- matrix(Mval, 1, lencol) 
	Mfree <- matrix(Mfree, 1, lencol)
	colnames(Mlab) <- colnames(Mval) <- colnames(Mfree) <- fcollab

	varzw <- 0
	if(!is.null(zwlab)) varzw <- var(as.vector(as.matrix(datawide[,zwlab]))) # Get the variance of the within-level covariate

	# These lines do nothing. Just prevent note from compiling packages.
	groupdiff <- NULL
	l1error <- NULL
	l2error <- NULL
	zbeffect <- NULL
	varzb <- NULL
	zweffect <- NULL
	
	constraint <- NULL # Make an appropriate constraint for each type of effect size and the existence of covariate
	if(estype == 0) { # Total variance to be standardized
		if(is.null(zwlab)) {
			if(is.null(zblab)) {
				constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error + l2error), name = "es")
			} else {
				constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error + l2error + (zbeffect^2 * varzb)), name = "es")
			}
		} else {
			if(is.null(zblab)) {
				constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error + (zweffect^2 * varzw) + l2error), name = "es")
			} else {
				constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error + (zweffect^2 * varzw) + l2error + (zbeffect^2 * varzb)), name = "es")
			}
		}
	} else if (estype == 1) { # Within-level variance to be standardized
		if(is.null(zwlab)) {
			constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error), name = "es")
		} else {
			constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error + (zweffect^2 * varzw)), name = "es")
		}	
	} else if (estype == 2) { # Between-level variance to be standardized
		if(is.null(zblab)) {
			constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l2error), name = "es")
		} else {
			constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l2error + (zbeffect^2 * varzb)), name = "es")
		}	
	} else {
		stop("'estype' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
	}
	onecov <- OpenMx::mxModel("effect size CRD", type="RAM",
	  OpenMx::mxData(datawide, type="raw"),
	  OpenMx::mxMatrix(type="Full", nrow=lencol, ncol=lencol, values=Aval, free=Afree, labels=Alab, name="A"),
	  OpenMx::mxMatrix(type="Symm", nrow=lencol, ncol=lencol, values=Sval, free=Sfree, labels=Slab, name="S"), 
	  OpenMx::mxMatrix(type="Full", nrow=lenrow, ncol=lencol, values=Fval, free=Ffree, labels=Flab, name="F"),
	  OpenMx::mxMatrix(type="Full", nrow=1, ncol=lencol, values=Mval, free=Mfree, labels=Mlab, name="M"), 
	  OpenMx::mxMatrix(type="Full", nrow=1, ncol=1, values=varzw, free=FALSE, labels="varzw", name="J"), 
	  OpenMx::mxRAMObjective("A","S","F","M", dimnames=fcollab), constraint, OpenMx::mxCI(c("es"), interval = conf.level)
	) 

	onecovfit <- OpenMx::mxRun(onecov, intervals=TRUE) # Run a Mx model with interval estimation with likelihood-based CI
	return(onecovfit@output$confidenceIntervals)
}


# Calculate a likelihood-based CI of ES
.likCIESCRDunequal <- function(datawide, ylab, xlab, zwlab=NULL, zblab=NULL, estype=1, iccy=0.25, es=0.5, totalvar=1, covariate=FALSE, iccz=0.25, r2within=0.5, r2between=0.5, totalvarz = 1, conf.level = 0.95) 
{
if(!requireNamespace("OpenMx", quietly = TRUE)) stop("The package 'OpenMx' is needed; please install the package and try again.")
  
	tau <- iccy * totalvar # Between-group variance
	sigma <- (1 - iccy) * totalvar # Within-group variance
	gamma1 <- NULL # Unstandardized conditions difference
	if(estype == 0) { # Create the unstandardized conditions difference based on the type of effect size
		gamma1 <- es * sqrt(totalvar)
	} else if (estype == 1) {
		gamma1 <- es * sqrt(sigma)
	} else if (estype == 2) {
		gamma1 <- es * sqrt(tau)
	} else {
		stop("'estype' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
}
	
if(covariate) 
	  {
		if(iccz == 0 && r2between != 0) stop("Because the covariate varies at the level 1 only, the r-square at level 2 must be 0.")
		if(iccz == 1 && r2within != 0) stop("Because the covariate varies at the level 2 only, the r-square at level 1 must be 0.")
		tauz <- totalvarz * iccz
		sigmaz <- totalvarz * (1 - iccz)
		gammazb <- sqrt(r2between * tau / tauz)
		gammazw <- sqrt(r2within * sigma / sigmaz)
		tau <- (1 - r2between) * tau # Update residual between-group variance
		sigma <- (1 - r2within) * sigma # Update residual within-group variance
}
	
	ntreat <- sum(sapply(datawide, function(x) sum(x[,xlab])))
	totaln <- sum(sapply(datawide, nrow))
	probx <- ntreat/totaln # Probability of the treatment condition
	
	FUNgroupsize <- function(dat, y, zw = NULL) {
		latentlab <- c("intcept", "slope") # Intercept = Between-group variation of dv; slope = The within-level effect of the covariate
		if(is.null(zw)) latentlab <- "intcept"

		# Get the dimensions of matrix in the RAM notation
		frowlab <- c(y, xlab, zblab) # Manifest variables labels
		fcollab <- c(frowlab, latentlab) # Manifest+Latent variables labels
		lenrow <- length(frowlab) 
		lencol <- length(fcollab)
		
		# Regression coefficients matrix
		Alab <- matrix(NA, lencol, lencol)
		Aval <- matrix(0, lencol, lencol)
		Afree <- matrix(FALSE, lencol, lencol)
		colnames(Alab) <- colnames(Aval) <- colnames(Afree) <- rownames(Alab) <- rownames(Aval) <- rownames(Afree) <- fcollab
		
		Alab["intcept", xlab] <- "groupdiff" # unstandardized conditions difference
		if(!is.null(zblab)) Alab["intcept", zblab] <- "zbeffect" # The between-level effect of the covariate
		if(!is.null(zw)) Alab[y, "slope"] <- paste("data.", zw, sep="") # Definition variables putting the within-level covariate values in the factor loadings from the "slope" factor
		
		Aval["intcept", xlab] <- gamma1 
		if(!is.null(zblab)) Aval["intcept", zblab] <- gammazb 
		Aval[y, latentlab] <- 1
		
		Afree["intcept", xlab] <- TRUE
		if(!is.null(zblab)) Afree["intcept", zblab] <- TRUE

		# Covariance matrix
		Slab <- matrix(NA, lencol, lencol)
		Sval <- matrix(0, lencol, lencol)
		Sfree <- matrix(FALSE, lencol, lencol)
		colnames(Slab) <- colnames(Sval) <- colnames(Sfree) <- rownames(Slab) <- rownames(Sval) <- rownames(Sfree) <- fcollab
		diag(Slab)[1:length(y)] <- "l1error" # Resiual variance in the within level
		diag(Sval)[1:length(y)] <- sigma 
		diag(Sfree)[1:length(y)] <- TRUE
		Slab["intcept", "intcept"] <- "l2error" # Residual variance in the between level
		Sval["intcept", "intcept"] <- tau 
		Sfree["intcept", "intcept"] <- TRUE
		Slab[xlab, xlab] <- "varx" # The variance of the grouping variables
		Sval[xlab, xlab] <- probx * (1 - probx) 
		Sfree[xlab, xlab] <- TRUE
		if(!is.null(zblab)) {
			Slab[c(xlab, zblab), c(xlab, zblab)] <- "covxzb" # The covariance between between-level covariate and grouping variable
			Slab[xlab, xlab] <- "varx" # The variance of the grouping variable
			Slab[zblab, zblab] <- "varzb" # The variance of the between-level covariate
			Sval[c(xlab, zblab), c(xlab, zblab)] <- 0
			Sval[xlab, xlab] <- probx * (1 - probx) 
			Sval[zblab, zblab] <- tauz 
			Sfree[c(xlab, zblab), c(xlab, zblab)] <- TRUE
		}

		# The selection matrix
		Fval <- cbind(diag(lenrow), matrix(0, lenrow, length(latentlab)))
		Flab <- matrix(NA, lenrow, lencol)
		Ffree <- matrix(FALSE, lenrow, lencol)
		colnames(Flab) <- colnames(Fval) <- colnames(Ffree) <- fcollab
		rownames(Flab) <- rownames(Fval) <- rownames(Ffree) <- frowlab

		# The intercept vectors
		Mlab <- c(rep(NA, length(y)), "meanX") # The proportion of treatment group
		Mval <- c(rep(0, length(y)), probx) 
		Mfree <- c(rep(FALSE, length(y)), TRUE)
		
		if(!is.null(zblab)) {
			Mlab <- c(Mlab, "meanzb") # The average of the between-level covariate
			Mval <- c(Mval, 0)
			Mfree <- c(Mfree, TRUE)	
		}
		Mlab <- c(Mlab, "meanctrl") # The average of the control condition
		Mval <- c(Mval, 0)
		Mfree <- c(Mfree, TRUE)
		if(!is.null(zw)) {
			Mlab <- c(Mlab, "zweffect") # The average of the within-level covariate
			Mval <- c(Mval, gammazw) 
			Mfree <- c(Mfree, TRUE)	
		} 	
		
		Mlab <- matrix(Mlab, 1, lencol) # Change a vector to an one-row matrix
		Mval <- matrix(Mval, 1, lencol) 
		Mfree <- matrix(Mfree, 1, lencol)
		colnames(Mlab) <- colnames(Mval) <- colnames(Mfree) <- fcollab

		onecov <- OpenMx::mxModel(paste0("group", length(y)), type="RAM",
		  OpenMx::mxData(dat, type="raw"),
			OpenMx::mxMatrix(type="Full", nrow=lencol, ncol=lencol, values=Aval, free=Afree, labels=Alab, name="A"),
			OpenMx::mxMatrix(type="Symm", nrow=lencol, ncol=lencol, values=Sval, free=Sfree, labels=Slab, name="S"), 
			OpenMx::mxMatrix(type="Full", nrow=lenrow, ncol=lencol, values=Fval, free=Ffree, labels=Flab, name="F"),
			OpenMx::mxMatrix(type="Full", nrow=1, ncol=lencol, values=Mval, free=Mfree, labels=Mlab, name="M"), 
			OpenMx::mxMatrix(type="Full", nrow=1, ncol=1, values=varzw, free=FALSE, labels="varzw", name="J"), 
			OpenMx::mxRAMObjective("A","S","F","M", dimnames=fcollab)
		)
		return(onecov)
	}
	varzw <- 0
	if(!is.null(zwlab)) varzw <- weighted.mean(do.call(c, mapply(function(x, y) var(as.vector(as.matrix(x[,y]))), x = datawide, y = zwlab,SIMPLIFY=FALSE)), as.numeric(names(datawide)))

	
	
	# These lines do nothing. Just prevent note from compiling packages.
	groupdiff <- NULL
	l1error <- NULL
	l2error <- NULL
	zbeffect <- NULL
	varzb <- NULL
	zweffect <- NULL
	
	constraint <- NULL # Make an appropriate constraint for each type of effect size and the existence of covariate
	if(estype == 0) { # Total variance to be standardized
		if(is.null(zwlab)) {
			if(is.null(zblab)) {
				constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error + l2error), name = "es")
			} else {
				constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error + l2error + (zbeffect^2 * varzb)), name = "es")
			}
		} else {
			if(is.null(zblab)) {
				constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error + (zweffect^2 * varzw) + l2error), name = "es")
			} else {
				constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error + (zweffect^2 * varzw) + l2error + (zbeffect^2 * varzb)), name = "es")
			}
		}
	} else if (estype == 1) { # Within-level variance to be standardized
		if(is.null(zwlab)) {
			constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error), name = "es")
		} else {
			constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l1error + (zweffect^2 * varzw)), name = "es")
		}	
	} else if (estype == 2) { # Between-level variance to be standardized
		if(is.null(zblab)) {
			constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l2error), name = "es")
		} else {
			constraint <- OpenMx::mxAlgebra(expression = groupdiff/sqrt(l2error + (zbeffect^2 * varzb)), name = "es")
		}	
	} else {
		stop("'estype' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
	}
	listModel <- NULL
	if(!is.null(zwlab)) {
		listModel <- mapply(FUNgroupsize, dat=datawide, y=ylab, zw=zwlab)
	} else {
		listModel <- mapply(FUNgroupsize, dat=datawide, y=ylab)
	}
	title <- "Effect Size CRD"
	algebra <- OpenMx::mxAlgebra("", name="allobjective")
	groupnames <- paste0("group", names(datawide))
	groupnames <- paste0(groupnames, ".objective")
	groupnames <- lapply(groupnames, as.name)
	algebra@formula <- as.call(c(list(as.name("sum")), groupnames))
	objective <- OpenMx::mxAlgebraObjective("allobjective")
	finalmodel <- OpenMx::mxModel(title, OpenMx::mxMatrix(type="Full", nrow=1, ncol=1, values=varzw, free=FALSE, labels="varzw", name="J"), unlist(listModel), constraint, algebra, objective, OpenMx::mxCI(c("es"), interval = conf.level))	
	finalmodelfit <- OpenMx::mxRun(finalmodel, intervals=TRUE)
	return(finalmodelfit@output$confidenceIntervals)
}

# Create data and find the width of the likelihood-based CI of ES
.runrepWidthESCRD <- function(seed, nclus, ntreatclus, nindiv, iccy, es, estype = 1, totalvar=1, covariate=FALSE, iccz=NULL, r2within=NULL, r2between=NULL, totalvarz = 1, conf.level = 0.95, diffsize = NULL) {
	set.seed(seed)
	datawide <- .createDataCRDWide(nclus=nclus, ntreatclus=ntreatclus, nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=covariate, iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = totalvarz, diffsize=diffsize) # Create data
	ylab <- NULL
	if(!is.null(diffsize)) {
		size <- as.numeric(names(datawide))
		ylab <- lapply(size, function(x) paste("y", 1:x, sep=""))
	} else {
		ylab <- paste("y", 1:nindiv, sep="")
	}

	xlab <- "x"
	zwlab <- NULL
	if(covariate && iccz != 1) {
		if(!is.null(diffsize)) {
			size <- as.numeric(names(datawide))
			zwlab <- lapply(size, function(x) paste("zw", 1:x, sep=""))
		} else {
			zwlab <- paste("zw", 1:nindiv, sep="")
		}
	}
	zblab <- NULL
	if(covariate && iccz != 0) zblab <- "zb"
	if(!is.null(diffsize)) {
		screencapture <- capture.output(
			result <- .likCIESCRDunequal(datawide=datawide, ylab=ylab, xlab=xlab, zwlab=zwlab, zblab=zblab, estype=estype, iccy=iccy, es=es, totalvar=totalvar, covariate=covariate, iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = totalvarz, conf.level = conf.level)) # Find the likelihood-based CI
	} else {
		screencapture <- capture.output(
			result <- .likCIESCRD(datawide=datawide, ylab=ylab, xlab=xlab, zwlab=zwlab, zblab=zblab, estype=estype, iccy=iccy, es=es, totalvar=totalvar, covariate=covariate, iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = totalvarz, conf.level = conf.level)
		) # Find the likelihood-based CI
	}
	return(result[2] - result[1])
}

# Simulate multiple datasets and find the average of the width of CI of ES (with or without the degree of assurance)
.findWidthCRDES <- function(nrep, nclus, ntreatclus, nindiv, iccy, es, estype = 1, totalvar=1, covariate=FALSE, iccz=NULL, r2within=NULL, r2between=NULL, totalvarz = 1, assurance=NULL, seed=123321, multicore=FALSE, numProc=NULL, conf.level=0.95, diffsize = NULL) {
	set.seed(seed)
	
	# Create list of seed number
	seedList <- as.list(sample(1:999999, nrep))
	Result.l <- NULL
    
	# Distribute the seed numbers and send to the runrepWidthESCRD function; Make the possibility of the multiple processors
    if (multicore) {
      if(!requireNamespace("parallel", quietly = TRUE)) stop("The package 'parallel' is needed; please install the package and try again.")
      sys <- .Platform$OS.type
      if (is.null(numProc)) 
        numProc <- parallel::detectCores()
      if (sys == "windows") {
        cl <- parallel::makeCluster(rep("localhost", numProc), type = "SOCK")                        
        Result.l <- parallel::clusterApplyLB(cl, seedList, .runrepWidthESCRD, nclus=nclus, ntreatclus=ntreatclus, nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=covariate, iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = totalvarz, conf.level = conf.level, diffsize=diffsize)
        parallel::stopCluster(cl)
      } else {
        Result.l <- parallel::mclapply(seedList, .runrepWidthESCRD, nclus=nclus, ntreatclus=ntreatclus, nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=covariate, iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = totalvarz, conf.level = conf.level, diffsize=diffsize)
      }
    } else {
      Result.l <- lapply(seedList, .runrepWidthESCRD, nclus=nclus, ntreatclus=ntreatclus, nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=covariate, iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = totalvarz, conf.level=conf.level, diffsize=diffsize)
    }
	result <- do.call(c, Result.l) # Change into a vector of width
	if(is.null(assurance)) { # Provide the average width or width of a specified degree of assurance
		return(mean(result, na.rm=TRUE))
	} else {
		return(quantile(result, assurance, na.rm=TRUE))
	}
}


# Find the number of clusters given the specified width of ES and the cluster size
.findNclusCRDES <- function(width, nindiv, es, estype = 1, iccy, prtreat, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, numProc=NULL, diffsize=NULL) {
	if(numpredictor > 0 & is.null(iccz)) iccz <- iccy
	if(numpredictor > 1) stop("Only one predictor is allowed.")
	totalvar <- 1
	if(estype == 0) {
		totalvar <- 1
	} else if (estype == 1) {
		totalvar <- 1/(1 - iccy)
	} else if (estype == 2) {
		totalvar <- 1/iccy
	} else {
		stop("'estype' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
	}
	startval <- .findNclusCRDDiff(width=width, nindiv=nindiv, prtreat=prtreat, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	startval <- as.numeric(startval)
	startwidth <- .findWidthCRDES(nrep, assurance=assurance, nclus=startval, ntreatclus=round(startval * prtreat), nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=as.logical(numpredictor), iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = 1, seed=seed, multicore=multicore, numProc=numProc, conf.level=conf.level, diffsize = diffsize)
	if(startwidth < width) {
		repeat {
			startval <- startval - 1
			if(round(startval * prtreat) == 1 | (startval - round(startval * prtreat)) == 1) return(c(startval + 1, startwidth))
			savedwidth <- startwidth
			startwidth <- .findWidthCRDES(nrep, assurance=assurance, nclus=startval, ntreatclus=round(startval * prtreat), nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=as.logical(numpredictor), iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = 1, seed=seed, multicore=multicore, numProc=numProc, conf.level=conf.level, diffsize = diffsize)
			if(startwidth > width) return(c(startval + 1, savedwidth))
		}
	} else if (startwidth > width) {
		repeat {
			startval <- startval + 1
			startwidth <- .findWidthCRDES(nrep, assurance=assurance, nclus=startval, ntreatclus=round(startval * prtreat), nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=as.logical(numpredictor), iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = 1, seed=seed, multicore=multicore, numProc=numProc, conf.level=conf.level, diffsize = diffsize)
			if(startwidth < width) return(c(startval, startwidth))
		}
	} else {
		return(c(startval, startwidth))
	}
}

# Find the cluster size given the specified width of ES and the number of clusters
.findNindivCRDES <- function(width, nclus, es, estype = 1, iccy, prtreat, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, numProc=NULL, diffsize=NULL) {
	if(numpredictor > 0 & is.null(iccz)) iccz <- iccy
	if(numpredictor > 1) stop("Only one predictor is allowed.")
	totalvar <- 1
	if(estype == 0) {
		totalvar <- 1
	} else if (estype == 1) {
		totalvar <- 1/(1 - iccy)
	} else if (estype == 2) {
		totalvar <- 1/iccy
	} else {
		stop("'estype' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
	}
	startval <- .findNindivCRDDiff(width=width, nclus=nclus, prtreat=prtreat, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	if(startval == "> 100000") stop("The starting number of individuals is > 100,000. With the specified number of clusters, it seems impossible to get the specified width.")
	startval <- as.numeric(startval)
	startwidth <- .findWidthCRDES(nrep, assurance=assurance, nclus=nclus, ntreatclus=round(nclus * prtreat), nindiv=startval, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=as.logical(numpredictor), iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = 1, seed=seed, multicore=multicore, numProc=numProc, conf.level=conf.level, diffsize = diffsize)
	if(startwidth < width) {
		repeat {
			startval <- startval - 1
			if(startval == 1) return(c(startval + 1, startwidth))
			savedwidth <- startwidth
			startwidth <- .findWidthCRDES(nrep, assurance=assurance, nclus=nclus, ntreatclus=round(nclus * prtreat), nindiv=startval, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=as.logical(numpredictor), iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = 1, seed=seed, multicore=multicore, numProc=numProc, conf.level=conf.level, diffsize = diffsize)
			if(startwidth > width) return(c(startval + 1, savedwidth))
		}
	} else if (startwidth > width) {
		repeat {
			startval <- startval + 1
			startwidth <- .findWidthCRDES(nrep, assurance=assurance, nclus=nclus, ntreatclus=round(nclus * prtreat), nindiv=startval, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=as.logical(numpredictor), iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = 1, seed=seed, multicore=multicore, numProc=numProc, conf.level=conf.level, diffsize = diffsize)
			if(startwidth < width) return(c(startval, startwidth))
		}
	} else {
		return(c(startval, startwidth))
	}
}

# Find the least expensive combination of the number of clusters and cluster size given the specified width of ES
.findMinCostCRDES <- function(width, cluscost=0, indivcost=1, es, estype = 1, iccy, prtreat, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, numProc=NULL, diffsize=NULL) {
	if(numpredictor > 0 & is.null(iccz)) iccz <- iccy
	if(numpredictor > 1) stop("Only one predictor is allowed.")
	totalvar <- 1
	if(estype == 0) {
		totalvar <- 1
	} else if (estype == 1) {
		totalvar <- 1/(1 - iccy)
	} else if (estype == 2) {
		totalvar <- 1/iccy
	} else {
		stop("'estype' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
	}
	startval <- .findMinCostCRDDiff(width=width, cluscost=cluscost, indivcost=indivcost, prtreat=prtreat, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)
	startval <- as.numeric(startval)
	startnindiv <- c(startval[2] - 1, startval[2], startval[2] + 1)
	result <- sapply(startnindiv, .findNclusCRDES, width=width, es=es, estype = estype, iccy = iccy, prtreat=prtreat, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance = assurance, conf.level = conf.level, nrep = nrep, iccz = iccz, seed = seed, multicore = multicore, numProc = numProc, diffsize = diffsize)
	resultnclus <- result[1,]
	resultwidth <- result[2,]
	startbudget <- mapply(.costCRD, nclus=resultnclus, nindiv=startnindiv, MoreArgs=list(cluscost=cluscost, indivcost=indivcost, diffsize = diffsize))
	if(which(startbudget == min(startbudget)) == 1) {
		repeat {
			startnindiv <- startnindiv - 1
			resultnclus[2:3] <- resultnclus[1:2]
			startbudget[2:3] <- startbudget[1:2]
			resultwidth[2:3] <- resultwidth[1:2]
			result <- .findNclusCRDES(width=width, nindiv=startnindiv[1], es=es, estype = estype, iccy = iccy, prtreat=prtreat, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance = assurance, conf.level = conf.level, nrep = nrep, iccz = iccz, seed = seed, multicore = multicore, numProc = numProc, diffsize = diffsize)
			resultnclus[1] <- result[1]
			resultwidth[1] <- result[2]
			startbudget[1] <- .costCRD(nclus=resultnclus[1], nindiv=startnindiv[1], cluscost=cluscost, indivcost=indivcost, diffsize = diffsize)
			if(which(startbudget == min(startbudget)) != 1) return(c(resultnclus[2], startnindiv[2], startbudget[2], resultwidth[2]))
		}
	} else if (which(startbudget == min(startbudget)) == 3) {
		repeat {
			startnindiv <- startnindiv + 1
			resultnclus[1:2] <- resultnclus[2:3]
			startbudget[1:2] <- startbudget[2:3]
			resultwidth[1:2] <- resultwidth[2:3]
			result <- .findNclusCRDES(width=width, nindiv=startnindiv[3], es=es, estype = estype, iccy = iccy, prtreat=prtreat, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance = assurance, conf.level = conf.level, nrep = nrep, iccz = iccz, seed = seed, multicore = multicore, numProc = numProc, diffsize = diffsize)
			resultnclus[3] <- result[1]
			resultwidth[3] <- result[2]
			startbudget[3] <- .costCRD(nclus=resultnclus[3], nindiv=startnindiv[3], cluscost=cluscost, indivcost=indivcost, diffsize = diffsize)
			if(which(startbudget == min(startbudget)) != 3) return(c(resultnclus[2], startnindiv[2], startbudget[2], resultwidth[2]))
		}	
	} else {
		return(c(resultnclus[2], startnindiv[2], startbudget[2], resultwidth[2]))
	}
}

# Find the least width of ES combination of the number of clusters and cluster size given the budget
.findMinWidthCRDES <- function(budget, cluscost=0, indivcost=1, es, estype = 1, iccy, prtreat, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, numProc=NULL, diffsize = NULL) {
	if(numpredictor > 0 & is.null(iccz)) iccz <- iccy
	if(numpredictor > 1) stop("Only one predictor is allowed.")
	totalvar <- 1
	if(estype == 0) {
		totalvar <- 1
	} else if (estype == 1) {
		totalvar <- 1/(1 - iccy)
	} else if (estype == 2) {
		totalvar <- 1/iccy
	} else {
		stop("'estype' can be 0 (total variance), 1 (level-1 variance), or 2 (level-2 variance) only.")
	}
	FUN <- function(nclus, nindiv) {
		.findWidthCRDES(nrep=nrep, assurance=assurance, nclus=nclus, ntreatclus=round(nclus * prtreat), nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=totalvar, covariate=as.logical(numpredictor), iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = 1, seed=seed, multicore=multicore, numProc=numProc, conf.level=conf.level, diffsize=diffsize)
	}
	startval <- .findMinWidthCRDDiff(budget=budget, cluscost=cluscost, indivcost=indivcost, prtreat=prtreat, totalvar=totalvar, iccy=iccy, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, diffsize = diffsize)

	startnclus <- c(startval[1] - 1, startval[1], startval[1] + 1)
	resultnindiv <- sapply(startnclus, .findNindivCRDBudget, budget=budget, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	startwidth <- mapply(FUN, nclus = startnclus, nindiv=resultnindiv)

	if(which(startwidth == min(startwidth)) == 1) {
		repeat {
			startnclus <- startnclus - 1
			resultnindiv[2:3] <- resultnindiv[1:2]
			startwidth[2:3] <- startwidth[1:2]
			resultnindiv[1] <- .findNindivCRDBudget(startnclus[1], budget=budget, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
			startwidth[1] <- FUN(nclus=startnclus[1], nindiv=resultnindiv[1])
			if(which(startwidth == min(startwidth)) != 1) return(c(startnclus[2], resultnindiv[2], startwidth[2]))
		}
	} else if (which(startwidth == min(startwidth)) == 3) {
		repeat {
			startnclus <- startnclus + 1
			resultnindiv[1:2] <- resultnindiv[2:3]
			startwidth[1:2] <- startwidth[2:3]
			resultnindiv[3] <- .findNindivCRDBudget(startnclus[3], budget=budget, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
			startwidth[3] <- FUN(nclus=startnclus[3], nindiv=resultnindiv[3])
			if(which(startwidth == min(startwidth)) != 3) return(c(startnclus[2], resultnindiv[2], startwidth[2]))
		}	
	} else {
		return(c(startnclus[2], resultnindiv[2], startwidth[2]))
	}
}


ss.aipe.crd.es.nclus.fixedwidth <- function(width, nindiv, es, estype = 1, iccy, prtreat, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, numProc=NULL, cluscost=NULL, indivcost=NULL, diffsize=NULL) {
	suppressWarnings(result <- .findNclusCRDES(width=width, nindiv=nindiv, es=es, estype = estype, iccy=iccy, prtreat=prtreat, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, nrep = nrep, iccz = iccz, seed = seed, multicore = multicore, numProc=numProc, diffsize=diffsize))
	calculatedCost <- NULL
	if(!is.null(cluscost) && !is.null(indivcost)) calculatedCost <- .costCRD(result[1], nindiv, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	.reportCRD(result[1], nindiv, result[2], cost=calculatedCost, es=TRUE, estype=estype, assurance=assurance, diffsize=diffsize) 
	invisible(result[1])
}


ss.aipe.crd.es.nindiv.fixedwidth <- function(width, nclus, es, estype = 1, iccy, prtreat, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, numProc=NULL, cluscost=NULL, indivcost=NULL, diffsize=NULL) {
	suppressWarnings(result <- .findNindivCRDES(width=width, nclus=nclus, es=es, estype = estype, iccy=iccy, prtreat=prtreat, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, nrep = nrep, iccz = iccz, seed = seed, multicore = multicore, numProc=numProc, diffsize=diffsize))
	calculatedCost <- NULL
	if(!is.null(cluscost) && !is.null(indivcost)) calculatedCost <- .costCRD(nclus, result[1], cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	.reportCRD(nclus, result[1], result[2], cost=calculatedCost, es=TRUE, estype=estype, assurance=assurance, diffsize=diffsize) 
	invisible(result[1])
}


ss.aipe.crd.es.nclus.fixedbudget <- function(budget, nindiv, cluscost, indivcost, nrep=NULL, prtreat=NULL, iccy=NULL, es=NULL, estype = 1, numpredictor = 0, iccz=NULL, r2within=NULL, r2between=NULL, assurance=NULL, seed=123321, multicore=FALSE, numProc=NULL, conf.level=0.95, diffsize=NULL) {
	nclus <- .findNclusCRDBudget(budget=budget, nindiv=nindiv, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	calculatedWidth <- NULL
	if(!is.null(nrep) && !is.null(prtreat) && !is.null(nindiv) && !is.null(iccy)) {
		suppressWarnings(calculatedWidth <- .findWidthCRDES(nrep=nrep, nclus=nclus, ntreatclus=round(nclus * prtreat), nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=1, covariate=as.logical(numpredictor), iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = 1, assurance=assurance, seed=seed, multicore=multicore, numProc=numProc, conf.level=conf.level, diffsize=diffsize))
	}
	calculatedCost <- .costCRD(nclus, nindiv, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	.reportCRD(nclus, nindiv, calculatedWidth, cost=calculatedCost, es=TRUE, estype=estype, assurance=assurance, diffsize=diffsize) 
	invisible(nclus)
}


ss.aipe.crd.es.nindiv.fixedbudget <- function(budget, nclus, cluscost, indivcost, nrep=NULL, prtreat=NULL, iccy=NULL, es=NULL, estype = 1, numpredictor = 0, iccz=NULL, r2within=NULL, r2between=NULL, assurance=NULL, seed=123321, multicore=FALSE, numProc=NULL, conf.level=0.95, diffsize=NULL) {
	nindiv <- .findNindivCRDBudget(budget=budget, nclus=nclus, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	calculatedWidth <- NULL
	if(!is.null(nrep) && !is.null(prtreat) && !is.null(nindiv) && !is.null(iccy)) {
		suppressWarnings(calculatedWidth <- .findWidthCRDES(nrep=nrep, nclus=nclus, ntreatclus=round(nclus * prtreat), nindiv=nindiv, iccy=iccy, es=es, estype = estype, totalvar=1, covariate=as.logical(numpredictor), iccz=iccz, r2within=r2within, r2between=r2between, totalvarz = 1, assurance=assurance, seed=seed, multicore=multicore, numProc=numProc, conf.level=conf.level, diffsize=diffsize))
	}
	calculatedCost <- .costCRD(nclus, nindiv, cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	.reportCRD(nclus, nindiv, calculatedWidth, cost=calculatedCost, es=TRUE, estype=estype, assurance=assurance, diffsize=diffsize) 
	invisible(nindiv)
}


ss.aipe.crd.es.both.fixedbudget <- function(budget, cluscost=0, indivcost=1, es, estype = 1, iccy, prtreat, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, numProc=NULL, diffsize=NULL) {
	suppressWarnings(result <- .findMinWidthCRDES(budget=budget, cluscost=cluscost, indivcost=indivcost, es=es, estype = estype, iccy=iccy, prtreat=prtreat, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, nrep = nrep, iccz = iccz, seed = seed, multicore = multicore, numProc=numProc, diffsize=diffsize))
	calculatedCost <- .costCRD(result[1], result[2], cluscost=cluscost, indivcost=indivcost, diffsize=diffsize)
	.reportCRD(result[1], result[2], result[3], cost=calculatedCost, es=TRUE, estype=estype, assurance=assurance, diffsize=diffsize) 
	invisible(result[1:2])
}


ss.aipe.crd.es.both.fixedwidth <- function(width, cluscost=0, indivcost=1, es, estype = 1, iccy, prtreat, r2between = 0, r2within = 0, numpredictor = 0, assurance=NULL, conf.level = 0.95, nrep = 1000, iccz = NULL, seed = 123321, multicore = FALSE, numProc=NULL, diffsize=NULL) {
	suppressWarnings(result <- .findMinCostCRDES(width=width, cluscost=cluscost, indivcost=indivcost, es=es, estype = estype, iccy=iccy, prtreat=prtreat, r2between = r2between, r2within = r2within, numpredictor = numpredictor, assurance=assurance, conf.level = conf.level, nrep = nrep, iccz = iccz, seed = seed, multicore = multicore, numProc=numProc, diffsize=diffsize))
	.reportCRD(result[1], result[2], result[4], cost=result[3], es=TRUE, estype=estype, assurance=assurance, diffsize=diffsize) 
	invisible(result[1:2])
}
