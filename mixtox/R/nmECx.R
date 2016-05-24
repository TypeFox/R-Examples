nmECx <- function(model, param, effv, minx){
	#calculate effect concentrations using associated inverse function
	if (missing(model) || missing (param) || missing(effv) || missing(minx)) stop('argument missing')
	if (is.vector(param)) param <- t(param)
	effv <- sort(effv)
	
	if(min(effv) < 0){
		effv_neg <- effv[effv < 0]
		effv_pos <- effv[effv >= 0]
		len_neg <- length(effv_neg)
		len_pos <- length(effv_pos)

		ecx <- matrix(0, length(model), (len_neg * 2 + len_pos))
		left_name <- paste('ECL', effv_neg * 100, sep = '')
		right_name <- paste('ECR', effv_neg * 100, sep = '')
		if(len_pos > 0) pos_name <- paste('EC', effv_pos * 100, sep = '') else pos_name <- c()
		colName <- c(left_name, right_name, pos_name)
		
	}else{
		len_pos <- length(effv)
		len_neg <- 0
		ecx <- matrix(0, length(model), len_pos)
		colName <- paste('EC', effv * 100, sep = '')
	}
	colnames(ecx) <- colName
	if (is.null(rownames(param)))  rownames(ecx) <- model else rownames(ecx) <- rownames(param)
	
	b <- 10000 # upper bound.
	eps <- 1e-10
	for(i in seq(model)){
		a <- minx[i]
		if(model[i] == 'Brain_Consens') f <- paste('1 - (1 + ', param[i, 1], '* xx) / (1 + exp(', param[i, 2], '*', param[i, 3],') * xx^',param[i, 2],')', sep = '')
		if(model[i] == 'BCV') f <- paste('1 -', param[i, 1], '* (1 +', param[i, 2], '* xx) / (1 + (1 + 2 *', param[i, 2], '*', param[i, 3],') * (xx /', param[i, 3],')^',param[i, 4],')', sep = '')
		if(model[i] == 'Cedergreen') f <- paste('1 - (1 +', param[i, 1], '* exp(-1 / (xx^',param[i, 2],'))) / (1 + exp(',param[i, 3], '* (log(xx) - log(', param[i, 4],')))', sep = '')
		if(model[i] == 'Beckon') f <- paste('(',param[i, 1], '+ (1 - (',param[i, 1],') / (1 + (',param[i, 2], '/ xx)^',param[i, 3],'))) / (1 + (xx /', param[i, 4],')^',param[i, 5],')', sep = '')
		if(model[i] == 'Biphasic') f <- paste(param[i, 1], '-', param[i, 1], '/ (1 + 10^((xx -', param[i, 2],') *', param[i, 3],')) + (1 -', param[i, 1],') / (1 + 10^((',param[i, 4], '- xx) *', param[i, 5],'))', sep = '')
		if(model[i] == 'Hill_six') f <- paste('(', param[i, 3], '/ (1 + (', param[i, 1], '/ xx)^', param[i, 2], ')) * (', param[i, 6], '/ (1 + (', param[i, 4], '/ xx)^', param[i, 5], '))')
		for(j in seq(len_pos + len_neg)){
			# the right side of minx
			fun_body <- paste(f, '-', effv[j], sep = '')
			fun = function(xx) eval(parse(text = fun_body))
			root <- uniroot(fun, c(a, b), tol = eps)$root
			ecx[i, (len_neg + j)] <- root
		}
	
		if(len_neg > 0){
			# the left side of minx
			for(k in seq(len_neg)){
				fun_body <- paste(f, '-', effv_neg[k], sep = '')
				fun = function(xx) eval(parse(text = fun_body))
				root <- uniroot(fun, c(eps, a), tol = eps)$root
				ecx[i, k] <- root
			}
		}
		
	}
	return(ecx)
}
