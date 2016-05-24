ecaPred <- function(effv, sgl, mix, pctMix){
	## calculating CA predicted effects at particular effect concentrations
	
	dichotomy <- function(fun, a, b, eps){
		# find the root using dichotomy method
		# xx means independent variable x
		
		expr <- parse(text = fun)
		execFun <- function(xx){}
		body(execFun) <- expr
		
		while (abs(b - a) > eps){
			ab2 <- (a + b) / 2
			flag <- sign(execFun(ab2) * execFun(a))
			if (flag == 0){
				ab2 <- (a + b) / 2
				break
			}
			else if(flag == -1)
				b <- ab2
			else if (flag == 1)
				a <- ab2
		}
		
		ab2 <- (a + b) / 2
		return(ab2)
	}

	eConAdd <- function(model, param, conc){
		# vector to matrix in case
		if (is.vector(conc)) conc <- t(conc)
		fac <- ncol(conc)
		lev <- nrow(conc)
		caFun <- as.character(rep(0, lev))
		
		for (i in seq(lev)){
			## CA equation construction
			for (j in seq(fac)){
				#if (model[j] == "Hill_two")
				#	caFun[i] <- paste(caFun[i], '+', conc[i, j], '/', param[j, 1],  '* xx / (', param[j, 2], '- xx)', sep = '')
				#else 
				if (model[j] == "Hill")
					caFun[i] <- paste(caFun[i], '+', conc[i, j], '/', param[j, 1], '/ ((1 / xx - 1)^(1 /', param[j, 2], '))', sep = '')
				else if (model[j] == "Weibull")
					caFun[i] <- paste(caFun[i], '+', conc[i, j], '/10^((log(-log(1 - xx))-', param[j, 1], ') /', param[j, 2], ')', sep = '')
				else if (model[j] == "Logit")
					caFun[i] <- paste(caFun[i], '+', conc[i, j], '/ 10^((log(xx / (1 - xx)) - ', param[j, 1], ') /', param[j, 2], ')', sep = '')
				else if (model[j] == "BCW")
					caFun[i] <- paste(caFun[i], '+', conc[i, j], '/(((', param[j, 3] / param[j, 2], ') * (log(-log(1 - xx)) -', param[j, 1], ') + 1)^(', 1 / param[j, 3], '))', sep = '')
				else if (model[j] == "BCL")
					caFun[i] <- paste(caFun[i], '+', conc[i, j], '/(((', param[j, 3] / param[j, 2], ') * (-log((1 - xx) / xx) -', param[j, 1], ') + 1)^(1/', param[j, 3], '))', sep = '')
				else if (model[j] == "GL")
					caFun[i] <- paste(caFun[i], '+', conc[i, j], '/(10^((-log((1 / xx)^(1/', param[j, 3], ')-1) - ', param[j, 1], ') / ', param[j, 2], '))', sep = '')
			}
		}
		
		a <- 0.0001 
		b <- 0.9999
		eps <- 1e-5
		eca <- rep(0, lev)
		
		for (i in seq(lev)){
			fun <- paste(caFun[i], '- 1', sep = '')
			eca[i] <- dichotomy(fun, a, b, eps)	
		}
		
		return(eca)
	}

	## source('ECx.R')
	## source('dichotomy.R')
	if (missing(effv) || missing(sgl) || missing(mix) || missing(pctMix))  stop('argument missing')
	
	model.sgl <- sgl$model
	param.sgl <- sgl$param
	model.mix <- mix$model
	param.mix <- mix$param
	len <- length(effv)
	eca <- matrix(0, len, length(model.mix))
	
	for (i in seq(len)){
		ecx.mix <- ECx(model.mix, param.mix, effv[i])
		ecx.vec <- as.vector(t(ecx.mix))
		ecx.conc <- ecx.vec * pctMix
		ecai <- eConAdd(model.sgl, param.sgl, ecx.conc)
		eca[i, ] <- ecai
	}
	
	rownames(eca) <- paste('caPred', '-E', effv * 100, sep = '')
	colnames(eca) <- rownames(pctMix)
	return(eca)
}
