caPred <- function(model, param, mixType = c("acr", "eecr", "udcr"), effv, effPoints){
	## concentration addition prediction
	##source('ECx.R')

	concAdd <- function(pctEcx, effPoints){
		## concentration addition
		ecPoints <- ECx(model, param, effPoints)
		ca <- 1 / (t(pctEcx) %*% (1 / ecPoints))
		return(ca)
	}

	if (missing(model) || missing(param) || missing(mixType) || missing(effv)) stop('argument missing')
	
	if (length(model) >= 2){
		## at these effect points the effect concentrations will be predicted
		if(missing(effPoints)){		
			effPoints <- c(.025, .03, .05, .10, .15, .20, .25, .30, .35, .40, .45, .47, .50, .52, .55, .60, .65, .70, .75, .80, .85, .90)
		}
		if (mixType == 'eecr'){
			## equal effect concentration ratio
			ecx <- ECx(model, param, effv)
			num <- nrow(ecx)
			mixEcx <- colSums(ecx)
			
			if (length(effv) > 1) pctEcx <- ecx / t(replicate(num, mixEcx)) else pctEcx <- ecx / mixEcx
			rownames(pctEcx) <- rownames(ecx)
			ca <- concAdd(pctEcx, effPoints)
			rowName <- paste('ca.EE', effv * 100, sep = '')
			rownames(ca) <- rowName
			designTable <- NULL
			
		}else if (mixType == 'acr'){
			## arbitrary concentration ratio
			if(length(model) != length(effv)) stop('no consistence')
			pctEcx <- t(t(effv / sum(effv)))
			ca <- concAdd(pctEcx, effPoints)
			rownames(ca) <- 'ca.acr'
			designTable <- NULL
			
		}else if(mixType == 'udcr'){
			## uniform design concentration ratio
			## source('unidTab.R')
			fac <- length(model)
			lev <- length(effv)
			tab <- unidTab(lev, fac)$T
			
			if (length(dim(tab)) == 3)
				## use the first uniform table if many
				uniTable <- tab[, , 1]
			
			if(length(dim(tab)) == 2)
				uniTable <- tab
			ecx <- ECx(model, param, effv)
			ecxMix <- matrix(0, fac, lev)
			
			## uniform mixture construction
			for (i in seq(fac)){
				for (j in seq(lev)){
					k <- uniTable[j, i]
					ecxMix[i, j] <- ecx[i, k]
				}
			}
			
			mixEcx <- colSums(ecxMix)
			pctEcx <- ecxMix / t(replicate(fac, mixEcx))
			ca <- concAdd(pctEcx, effPoints)
			rowName <- paste('ca.U', seq(lev), sep = '')
			rownames(ca) <- rowName
			rownames(pctEcx) <- rownames(ecx)
			colnames(pctEcx) <- rowName
			designTable <- uniTable
		}

	}else {
		stop('needs more than one component')
	}
	list(ca = ca, e = effPoints, pct = t(pctEcx), unitab = designTable)
}