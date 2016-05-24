CEx <- function(model, param, conc){
	# calculate response based on concentration
	if (missing(model) || missing (param) || missing(conc)) stop('argument missing')
	#if (missing(conc)) conc = 0.00005
	if (is.vector(param)) param <- t(param)
	
	effv <- matrix(0, length(model), length(conc))
	
	for (i in seq(model)){
		fun <- model[i]
		p <- param[i, ]
		
		for (j in seq(conc)){
			if (fun == 'Hill')
				ev <- 1 / (1 + (p[1] / conc[j])^p[2])
			else if (fun == 'Hill_two')
				ev <- p[2] * conc[j] / (p[1] + conc[j])
			else if (fun == 'Hill_three')
				ev <- p[3] /(1 + (p[1] / conc[j])^p[2])
			else if (fun == 'Hill_four')
				ev <- p[4] + (p[3] - p[4]) / (1 + (p[1] / conc[j])^p[2])
			else if(fun == 'Weibull')
				ev <- 1 - exp(-exp(p[1] + p[2] * log10(conc[j])))
			else if(fun == 'Weibull_three')
				ev <- p[3] * (1 - exp(-exp(p[1] + p[2] * log10(conc[j]))))
			else if(fun == 'Weibull_four')
				ev <- p[3] + (p[4] - p[3]) * exp(-exp(p[1] + p[2] * log10(conc[j])))
			else if (fun == "Logit")
				ev <- 1 / (1 + exp(-p[1] - p[2] * log10(conc[j])))
			else if(fun == 'Logit_three')
				ev <- p[3] / (1 + exp((-p[1]) - p[2] * log10(conc[j])))
			else if(fun == 'Logit_four')
				ev <- p[4] + (p[3] - p[4]) / (1 + exp((-p[1]) - p[2] * log10(conc[j])))
			else if (fun == "BCW")
				ev <- 1 - exp(-exp(p[1] + p[2] * ((conc[j]^p[3] - 1) / p[3])))
			else if (fun == "BCL")
				ev <- 1 / (1 + exp(-p[1] - p[2]((conc[j]^p[3] - 1) / p[3])))
			else if (fun == "GL")
				ev <- 1 / (1 + exp(-p[1] - p[2] * log10(conc[j])))^p[3]
			else if (fun == "Brain_Consens") 
				ev <- 1 - (1 + p[1] * conc[j]) / (1 + exp(p[2] * p[3]) * conc[j]^p[2])
			else if(fun == "BCV") 
				ev <- 1 - p[1] * (1 + p[2] * conc[j]) / (1 + (1 + 2 * p[2] * p[3]) * (conc[j] / p[3])^p[4])
			else if(fun == "Cedergreen") 
				ev <- 1 - (1 + p[1] * exp(-1 / (conc[j]^p[2]))) / (1 + exp(p[3] * (log(conc[j]) - log(p[4]))))
			else if(fun == "Beckon") 
				ev <- (p[1] + (1 - (p[1]) / (1 + (p[2] / conc[j])^p[3]))) / (1 + (conc[j] / p[4])^p[5])
			else if(fun == "Biphasic") 
				ev <- p[1] - p[1] / (1 + 10^((conc[j] - p[2]) * p[3])) + (1 - p[1]) / (1 + 10^((p[4] - conc[j]) * p[5]))
			else if(fun == 'Hill_six')
				ev <- (p[3] / (1 + (p[1] / conc[j])^p[2])) * (p[6] / (1 + (p[4] / conc[j])^p[5]))
			effv[i, j] <- ev
		}
	}
	
	colName <- paste('Rspn_@_', conc, sep = '')
	colnames(effv) <- colName
	if(is.null(rownames(param))) rownames(effv) <- model else rownames(effv) <- rownames(param)
	return(effv)
}
