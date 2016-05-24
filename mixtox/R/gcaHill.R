gcaHill <- function(model, param, mixType = c("acr", "eecr", "udcr"), effv, refEffv = c(0.10, 0.50)){
	## generalized concentration addition based only on Hill_two function
	gconcAdd <- function(model, param, pctEcx, refEffv){
	# concentration addition
		
		pointNum <- 22
		#dilution = 20
		refEcx <- ECx(model, param, refEffv)
		refEcx_new <- refEcx[refEcx > 0]
		refMin <- min(refEcx_new)
		refMax <- max(refEcx_new)
		#conc <- 10^(seq(log10(refMin / dilution), log10(refMax * dilution), length.out = pointNum))
		conc <- 10^(seq(log10(refMin / 2), log10(refMax * 4), length.out = pointNum))
		fac <- nrow(pctEcx)
		lev <- ncol(pctEcx)
		emix <- matrix(0, lev, pointNum)
		
		for(i in seq(lev)){
			for(j in seq(pointNum)){
				emix[i, j] <- sum(param[, 2] * conc[j] * pctEcx / param[, 1]) / (1 + sum(conc[j] * pctEcx / param[, 1]))
			}
		}
		list(x = conc, y = emix)
	}
	
	if (missing(model) || missing(param) || missing(mixType) || missing(effv)) stop('argument missing')
	
	if (length(model) >= 2){
	
		if (mixType == 'eecr'){
			## equal effect concentration ratio
			ecx <- ECx(model, param, effv)
			num <- nrow(ecx)
			mixEcx <- colSums(ecx)
			if (length(effv) > 1) pctEcx <- ecx / t(replicate(num, mixEcx)) else pctEcx <- ecx / mixEcx
			rownames(pctEcx) <- rownames(ecx)
			gca <- gconcAdd(model, param, pctEcx, refEffv)
			rowName <- paste('gca.EE', effv, sep = '')
			rownames(gca$y) <- rowName
			designTable <- NULL
			
		}else if (mixType == 'acr'){
			## arbitrary concentration ratio
			if(length(model) != length(effv)) stop('no consistence')
			pctEcx <- t(t(effv / sum(effv)))
			gca <- gconcAdd(model, param, pctEcx, refEffv)
			rownames(gca$y) <- 'ca.acr'
			designTable <- NULL
			
		}else if(mixType == 'udcr'){
			# uniform design
			## source('unidTab.R')
			fac <- length(model)
			lev <- length(effv)
			tab <- unidTab(lev, fac)$T
			
			if (length(dim(tab)) == 3)
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
			gca <- gconcAdd(model, param, pctEcx, refEffv)
			rowName <- paste('ca.U', seq(lev), sep = '')
			rownames(gca$y) <- rowName
			rownames(pctEcx) <- rownames(ecx)
			colnames(pctEcx) <- rowName
			designTable <- uniTable
		}

	}else {
		stop('needs more than one component')
	}
	list(x = gca$x, e = gca$y, pct = t(pctEcx), unitab = designTable)
}
