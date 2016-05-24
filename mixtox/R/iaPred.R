iaPred <- function(model, param, mixType = c("acr", "eecr", "udcr"), effv, effPoints, lb = 1e-9, ub = 6){
	## independent action prediction
	## source('ECx.R')
	
	indAct <- function(model, param, pctEcx, effPoints){
		# independent action
		ecPoints <- ECx(model, param, effPoints)
		fac <- nrow(pctEcx)
		lev <- ncol(pctEcx)
		iaFun <- as.character(rep(1, lev))
		
		for (i in seq(lev)){
			# IA equation construction
			# xx means x elsewhere
			for (j in seq(fac)){
				#if (model[j] == 'Hill_two')
				#	iaFun[i] <- paste(iaFun[i], '*', '(1 - (', param[j, 2], '* xx / (', param[j, 1], '+ xx)))', sep = '')
				#else 
				if (model[j] == 'Hill')
					iaFun[i] <- paste(iaFun[i], '*', '(1 - ( 1 / (1 + (', param[j, 1], '/', pctEcx[j, i], '* xx)^', param[j, 2], ')))', sep = '')
				else if (model[j] == "Weibull")
					iaFun[i] <- paste(iaFun[i], '*', '(1 - (1 - exp(-exp(', param[j, 1], '+', param[j, 2], '* log10(', pctEcx[j, i], '* xx)))))', sep = '')
				else if (model[j] == "Logit")
					iaFun[i] <- paste(iaFun[i], '*', '(1 - (1 / (1 + exp(-', param[j, 1], '-', param[j, 2], '* log10(', pctEcx[j, i], '* xx)))))', sep = '')
				else if (model[j] == "BCW")
					iaFun[i] <- paste(iaFun[i], '*', '(1 - (1 - exp(-exp(', param[j, 1], '+', param[j, 2], '* (((', pctEcx[j, i], '* xx)^', param[j, 3], '- 1) /', param[j, 3], ')))))', sep = '')
				else if (model[j] == "BCL")
					iaFun[i] <- paste(iaFun[i], '*', '(1 - ((1 + exp(-', param[j, 1], '-', param[j, 2], '* (((', pctEcx[j, i], '* xx)^', param[j, 3], '- 1) /', param[j, 3], ')))^(-1)))', sep = '')
				else if (model[j] == "GL")
					iaFun[i] <- paste(iaFun[i], '*', '(1 - (1 / (1 + exp(-', param[j, 1], '-', param[j, 2], '* log10(', pctEcx[j, i], '* xx)))^', param[j ,3], '))', sep = '')
			}
		}
		
		a <- lb
		b <- ub
		eps <- 1e-10		
		root <- matrix(0, lev, ncol(ecPoints))
		
		for (i in seq(lev)){
			fia <-  iaFun[i]
			for (k in seq(ncol(ecPoints))){
				value <- 1 - effPoints[k]
				fun <- paste(value, '-',  fia, sep = '')
				f = function(xx) eval(parse(text = fun))
				root[i, k] <- uniroot(f, c(a, b), tol = eps)$root
			}
		}
		
		colName <- paste('EC', effPoints * 100, sep = '')
		colnames(root) <- colName
		return(root)
	}
	
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
			ia <- indAct(model, param, pctEcx, effPoints)
			rowName <- paste('ia.EE', effv * 100, sep = '')
			rownames(ia) <- rowName
			designTable <- NULL
		}else if (mixType == 'acr'){
			## arbitrary concentration ratio
			if(length(model) != length(effv)) stop('no consistence')
			pctEcx <- t(t(effv / sum(effv)))
			ia <- indAct(model, param, pctEcx, effPoints)
			rownames(ia) <- 'ia.acr'
			designTable <- NULL
		}else if(mixType == 'udcr'){
			## uniform design concentration ratio
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
			ia <- indAct(model, param, pctEcx, effPoints)
			rowName <- paste('ia.U', seq(lev), sep = '')
			rownames(ia) <- rowName
			rownames(pctEcx) <- rownames(ecx)
			colnames(pctEcx) <- rowName
			designTable <- uniTable
		}

	}else {
		stop('needs more than one component')
	}
	
	list(ia = ia, e = effPoints, pct = t(pctEcx), unitab = designTable)
}