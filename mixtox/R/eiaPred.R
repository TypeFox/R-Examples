eiaPred <- function(effv, sgl, mix, pctMix) {
	## calculating IA predicted effects at particular effect concentrations
	
	eIndAct <- function(model, param, conc){
		# vector to matrix in case
		if (is.vector(conc)) conc <- t(conc)
		fac <- ncol(conc)
		lev <- nrow(conc)
		ia <- rep(1, lev)
		
		for (i in seq(lev)){
			## IA equation construction
			for (j in seq(fac)){
				#if (model[j] == "Hill_two")
				#	ia[i] <- ia[i] * (1 - (param[j, 2] * conc[i, j] / (param[j, 1] + conc[i, j])))
				#else 
				if (model[j] == "Hill")
					ia[i] <- ia[i] * (1 - 1 / (1 + (param[j, 1] / conc[i, j])^param[j, 2]))
				else if (model[j] == "Weibull")
					ia[i] <- ia[i] * (1 - (1 - exp(-exp(param[j, 1] + param[j, 2] * log10(conc[i, j])))))
				else if (model[j] == "Logit")
					ia[i] <- ia[i] * (1 - (1 / (1 + exp(-param[j, 1] - param[j, 2] * log10(conc[i, j])))))
				else if (model[j] == "BCW")
					ia[i] <- ia[i] * (1 - (1 - exp(-exp(param[j, 1] + param[j, 2] * ((conc[i, j]^param[j, 3] - 1) / param[j, 3])))))
				else if (model[j] == "BCL")
					ia[i] <- ia[i] * (1 - ((1 + exp(-param[j, 1] - param[j, 2] * ((conc[i, j]^param[j, 3] - 1) / param[j, 3])))^(-1)))
				else if (model[j] == "GL")
					ia[i] <- ia[i] * (1 - (1 / (1 + exp(-param[j, 1] - param[j, 2] * log10(conc[i, j])))^param[j, 3]))
			}
		}
		ia <- (1 - ia)
		return(ia)
	}
	
	## source('ECx.R')
	if (missing(effv) || missing(sgl) || missing(mix) || missing(pctMix))  stop('argument missing')
	
	model.sgl <- sgl$model
	param.sgl <- sgl$param
	model.mix <- mix$model
	param.mix <- mix$param
	len <- length(effv)
	eia <- matrix(0, len, length(model.mix))
	
	for (i in seq(len)){
		ecx.mix <- ECx(model.mix, param.mix, effv[i])
		ecx.vec <- as.vector(t(ecx.mix))
		ecx.conc <- ecx.vec * pctMix
		eiai <- eIndAct(model.sgl, param.sgl, ecx.conc)
		eia[i, ] <- eiai
	}
	
	rownames(eia) <- paste('iaPred', '-E', effv * 100, sep = '')
	colnames(eia) <- rownames(pctMix)
	return(eia)
}