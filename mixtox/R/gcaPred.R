gcaPred <- function(model, param, mixType = c("acr", "eecr", "udcr"), effv, refEffv = c(0.1, 0.50, 0.80), lb = 1E-7, ub = 0.9){
	# generalized concentration addition prediction
	#
	## source('ECx.R')
	
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


	gConcAdd <- function(model, param, pctEcx, refEffv, lb, ub){
		# generalized concentration addition
		#refEffv <- c(0.05, 0.50, 0.90)
		pointNum <- 22
		#refEffv <- 0.5
		#dilution <- 20

		refEcx <- ECx(model, param, refEffv)
		refEcx_new <- refEcx[refEcx > 0]
		refMin <- min(refEcx_new)
		refMax <- max(refEcx_new)
		conc <- 10^(seq(log10(refMin), log10(refMax), length.out = pointNum))
		#conc <- 10^(seq(log10(refMin / dilution), log10(refMax * dilution), length.out = pointNum))
		fac <- nrow(pctEcx)
		lev <- ncol(pctEcx)
		root <- matrix(0, lev, pointNum)
		
		for (i in seq(lev)){
			for(j in seq(pointNum)){
				fun <- as.character(1)
				for (k in seq(fac)){
				
					if (model[k] == 'Hill_two')
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '/ (', param[k, 1], '* xx / (', param[k, 2], '- xx))', sep = '')
					else if (model[k] == "Hill_three")
						# inv(y) = 1/Gamma*(1+(Alpha/x)^Beta)
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '/ (1 /', param[k, 3], '* (1 + (', param[k, 1], '/ xx)^', param[k, 2], '))', sep  = '')
					else if (model[k] == "Hill_four")
						# inv(y) = 1/(Delta+(Gamma-Delta)/(1+(Alpha/x)^Beta))
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '* (', param[k, 4], '+ (', param[k, 3], '-', param[k, 4], ') / (1 + (', param[k, 1], '/ xx)^', param[k, 2], '))', sep = '')
					else if (model[k] == "Weibull_three")
						#inv(y) = 1/Gamma/(1-exp(-exp(Alpha+Beta*log(x)/log(10))))
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '*', param[k, 3], '* (1 - exp(-exp(', param[k, 1], '+', param[k, 2], '*log(xx)/log(10))))', sep = '')
					else if (model[k] == "Weibull_four")
						#inv(y) = 1/(Gamma+(Delta-Gamma)*exp(-exp(Alpha+Beta*log(x)/log(10))))
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '*(', param[k, 3], '+(', param[k, 4], '-', param[k, 3], ')*exp(-exp(', param[k, 1], '+', param[k, 2], '* log(xx) / log(10))))', sep = '')
					else if (model[k] == "Logit_three")
						# inv(y) = 1/Gamma*(1+exp(-Alpha-Beta*log(x)/log(10)))
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '*', param[k, 3], '/ (1 + exp(-', param[k, 1], '-', param[k, 2], '* log(xx) / log(10)))', sep = '')
					else if (model[k] == "Logit_four")
						# inv(y) = 1/(Delta+(Gamma-Delta)/(1+exp(-Alpha-Beta*log(x)/log(10))))
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '*(', param[k, 4], '+(', param[k, 3], '-', param[k, 4], ')/(1 + exp(-', param[k, 1], '-', param[k, 2], '*log(xx) / log(10))))', sep = '')
					else if (model[k] == "Hill")
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '/', param[k, 1], '/ (1 / xx - 1)^(1 /', param[k, 2], ')', sep = '')
					else if (model[k] == "Weibull")
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '/ (exp(-(-log(log(-1 / (-1 + xx))) +', param[k, 1], ') * log(10) /', param[k, 2], '))', sep = '')
					else if (model[k] == "Logit")
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '/ (exp(-log(10) * (', param[k, 1], '+ log(-(-1 + xx) / (xx))) / ', param[k, 2], '))', sep = '')
					else if (model[k] == "BCW")
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '/(exp(log(-(', param[k, 1], '*', param[k, 3], '-', param[k, 2], ' - log(-log(1 - xx)) *', param[k, 3], ') /', param[k, 2], ') /', param[k, 3], '))', sep = '')
					else if (model[k] == "BCL")
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '/(exp(log(-(', param[k, 1], '*', param[k, 3], '-', param[k, 2], '+ log(-(-1 + xx) / (xx)) *', param[k, 3], ') /', param[k, 2], ') /', param[k, 3], '))', sep = '')
					else if (model[k] == "GL")
						fun <- paste(fun, '-', pctEcx[k, i] * conc[j], '/ (exp(-log(10) * (', param[k, 1], '+ log(exp(-log(xx) /', param[k, 3], ') - 1)) /', param[k, 2], '))', sep = '')
				}
				
				a <- lb
				b <- ub
				eps <- 1e-10
				#f = function(xx) eval(parse(text = fun))
				#root[i, j] <- uniroot(f, c(a, b), tol = eps)$root
				root[i, j] <- dichotomy(fun, a, b, eps)
			}
		}
		list(x = conc, y = root)
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
			gca <- gConcAdd(model, param, pctEcx, refEffv, lb, ub)
			rownames(gca$y) <- paste('gca.EE', effv, sep = '')
			designTable <- NULL
			
		}else if (mixType == 'acr'){
			## arbitrary concentration ratio
			
			if(length(model) != length(effv)) stop('no consistence')
			
			pctEcx <- t(t(effv / sum(effv)))
			gca <- gConcAdd(model, param, pctEcx, refEffv, lb, ub)
			rownames(gca$y) <- 'gca.acr'
			designTable <- NULL
			
		}else if(mixType == 'udcr'){
			## uniform design concentration ratio
			
			## source('unidTab.R')
			fac <- length(model)
			lev <- length(effv)
			tab <- unidTab(lev, fac)$T
			
			if (length(dim(tab)) == 3) uniTable <- tab[, , 1]
			
			if(length(dim(tab)) == 2) uniTable <- tab
			
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

			gca <- gConcAdd(model, param, pctEcx, refEffv, lb, ub)
			rowName <- paste('gca.U', seq(lev), sep = '')
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
