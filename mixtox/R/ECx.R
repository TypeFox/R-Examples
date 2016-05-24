ECx <- function(model, param, effv){
	#calculate effect concentrations using associated inverse function
	if (missing(model) || missing (param)) stop('argument missing')
	if (missing(effv)) effv = 0.5
	if (is.vector(param)) param <- t(param)
	
	effv <- sort(effv)
	ecx <- matrix(0, length(model), length(effv))

	for (i in seq(model)){
		fun <- model[i]
		p <- param[i, ]
		if (fun == 'Hill')
			ec <- p[1] / ((1 / effv - 1)^(1 / p[2]))
		else if (fun == 'Hill_two')
			ec <- p[1] * effv / (p[2] - effv)
		else if (fun == 'Hill_three')
			ec <- p[1] / ((p[3] / effv - 1)^(1 / p[2]))
		else if(fun == 'Hill_four')
			ec <- p[1] / (((p[3] - p[4]) / (effv - p[4]) - 1)^(1 / p[2]))
		else if(fun == 'Weibull')
			ec <- exp(-(-log(log(-1 / (-1 + effv))) + p[1]) * log(10) / p[2])
		else if(fun == 'Weibull_three')
			ec <- exp(-(-log(log(p[3] / (p[3] - effv))) + p[1]) * log(10) / p[2])
		else if(fun == 'Weibull_four')
			ec <- exp((log(log((-p[4] + p[3]) / (p[3] - effv))) - p[1]) * log(10) / p[2])
		else if (fun == "Logit")
			ec <- exp(-log(10) * (p[1] + log(-(-1 + effv) / (effv))) / p[2])
		else if(fun == 'Logit_three')
			ec <- exp(-log(10) * (p[1] + log((p[3] - effv) / effv)) / p[2])
		else if(fun == 'Logit_four')
			ec <- exp(-log(10) * (p[1] + log(-(p[3] - effv) / (p[4] - effv))) / p[2])
		else if (fun == "BCW")
			ec <- exp(log(-(p[1] * p[3] - p[2] - log(-log(1 - effv)) * p[3]) / p[2]) / p[3])
		else if (fun == "BCL")
			ec <- exp(log(-(p[1] * p[3] - p[2] + log(-(-1 + effv) / effv) * p[3]) / p[2]) / p[3])
		else if (fun == "GL")
			ec <- exp(-log(10) * (p[1] + log(exp(-log(effv) / p[3]) - 1)) / p[2])

		ecx[i, ] <- ec
	}
	
	colName <- paste('EC', effv * 100, sep = '')
	colnames(ecx) <- colName
	if (is.null(rownames(param)))  rownames(ecx) <- model else rownames(ecx) <- rownames(param)
	return(ecx)
}
