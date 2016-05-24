ds_test <- function(y, x, ..., type = c("ds", "eqp"), lambda = 1, alpha = 1, rounds = 0)
{
	if(length(y[is.na(y)]) > 0){
		stop("'y' should not contain NA.")
	}
	if(length(x[is.na(x)]) > 0){
		stop("'x' should not contain NA.")
	}
	type <- match.arg(type)
	DNAME <- deparse(substitute(y))
	if(lambda < 0){
		stop("Improper value of lambda, it should be larger than 0.")
	}else if(lambda == 0){
		warning("lambda equals zero. Be aware of that this means no penalty for number of slice.")
	}
	rounds <- floor(rounds)
	pvalue <- NULL
	if(is.integer(x)){
		if(length(x) != length(y)){
			stop("'x' and 'y' lengths differ.")
		}
		hist <- table(x)
		if(length(hist) < 2){
			stop("All the values of x are the same.")
		}
		if(max(hist) < 5){
			stop("Not enough data")
		}
		DNAME <- paste(DNAME, "and", deparse(substitute(x)))
		METHOD <- "K-sample test via dynamic slicing"
		xdim <- max(x) + 1
		x <- x[order(y)]
		n <- length(x)
		if(type == "ds"){
			dsres <- .Call('dslice_dslice_k', x, xdim, lambda, PACKAGE = 'dslice')
			if(rounds > 1){
				if(dsres$dsval > 1e-6){
					nullval <- vector(length = n, mode = "numeric")
					for(i in 1:rounds){
						nullval[i] <- .Call('dslice_ds_k', sample(x), xdim, lambda, PACKAGE = 'dslice')
					}
					pvalue <- length(which(nullval > dsres$dsval)) / rounds
				}else{
					pvalue <- 1
				}
			}
		}
		if(type == "eqp"){
			METHOD <- paste(METHOD, "with O(sqrt{n}) resolution", sep = " ")
			dsres <- .Call('dslice_dslice_eqp_k', x, xdim, lambda, PACKAGE = 'dslice')
			if(rounds > 1){
				if(dsres$dsval > 1e-6){
					nullval <- vector(length = n, mode = "numeric")
  					for(i in 1:rounds){
  						nullval[i] <- .Call('dslice_ds_eqp_k', sample(x), xdim, lambda, PACKAGE = 'dslice')
					}
					pvalue <- length(which(nullval > dsres$dsval)) / rounds
				}else{
					pvalue <- 1
				}
			}
		}
		STATISTIC <- dsres$dsval
		names(STATISTIC) <- "DS"
		ALTER <- "The distribution of Y given X = j (j = 1, ..., K) are not the same"
		RVAL <- list(statistic = STATISTIC, p.value = pvalue, alternative = ALTER, 
		             method = METHOD, data.name = DNAME, slices = dsres$slices)
	}else{
		if(is.character(x)){
			x <- get(x, mode = "function", envir = parent.frame())
		}
		if(!is.function(x)){
			stop("'x' must be numeric or a function or a string naming a valid function.")
		}
		METHOD <- "One-sample test via dynamic slicing"
		y <- sort(unique(y))
		n <- length(y)
		q <- x(y, ...)
		if(type == "ds"){
			if(alpha < 1){
				stop("Improper value of alpha, it should not less than 1.")
			}
			dsval <- .Call('dslice_ds_1', q, lambda, alpha, PACKAGE = 'dslice')
			if(rounds > 1){
				if(dsval > 1e-6){
					nullval <- vector(length = n, mode = "numeric")
					for(i in 1:rounds){
						nullval[i] <- .Call('dslice_ds_1', sort(runif(n, 0, 1)), lambda, alpha, PACKAGE = 'dslice')
					}
					pvalue <- length(which(nullval > dsval)) / rounds
				}else{
					pvalue <- 1
				}
			}
		}
		if(type == "eqp"){
			METHOD <- paste(METHOD, "with O(n) resolution", sep = " ")
			dsval <- .Call('dslice_ds_eqp_1', q, lambda, PACKAGE = 'dslice')
			if(rounds > 1){
				if(dsval > 1e-6){
					nullval <- vector(length = n, mode = "numeric")
					for(i in 1:rounds){
						nullval[i] <- .Call('dslice_ds_eqp_1', sort(runif(n, 0, 1)), lambda, PACKAGE = 'dslice')
					}
					pvalue <- length(which(nullval > dsval)) / rounds
				}else{
					pvalue <- 1
				}
			}
		}
		STATISTIC <- dsval
		names(STATISTIC) <- "DS"
		ALTER <- "Data are not drawn from null distribution."
		RVAL <- list(statistic = STATISTIC, p.value = pvalue, alternative = ALTER, 
		             method = METHOD, data.name = DNAME)
	}
	class(RVAL) <- "htest"
	return(RVAL)
}
