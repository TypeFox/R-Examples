budget.ie <- function(bdgt, method = "leave out", lo = 2, it = 100){
	METHODS <- c("leave out", "bootstrap", "sorted bootstrap", "constrained bootstrap", "jackknife", "jack-validate")
	method <- pmatch(method, METHODS)
	if (is.na(method)){
		stop("invalid method")
		}
	# define function
	do.budget <- function(z){
		m.sub <- models[z]
		bv <- update(bdgt, models=m.sub, set.back=NULL, return.models=FALSE)
		b <- sum(bv[,1], na.rm=TRUE)
		cat(paste(substitute(z)[[3]],":",sep=""))
		return(b)
	}
	# prepare
	models <- bdgt$models
	ind.year <- seq(length(models))
	# do
	if(method == 1){
		ind <- lapply(seq(it), function(z) sort(sample(ind.year, length(ind.year)-lo)))
		res <- sapply(ind, do.budget)
	}
	else if (method == 2){
		ind <- lapply(seq(it), function(z) sample(seq(length(ind.year)), replace=TRUE))
		res <- sapply(ind, do.budget)
	}
	else if (method == 3){
		ind <- lapply(seq(it), function(z) sort(sample(seq(length(ind.year)), replace=TRUE)))
		res <- sapply(ind, do.budget)
	}
	else if (method == 4){
		ind <- lapply(seq(it), function(z) unique(sort(sample(seq(length(ind.year)), replace=TRUE))))
		res <- sapply(ind, do.budget)
	}	
	else if (method == 5){
		ind <- combs(seq(length(ind.year)), length(ind.year)-1)
		res <- apply(ind, 1, do.budget)
	}	
	else if (method == 6){
		ind <- lapply(c(1:lo), function(z) combs(seq(length(ind.year)), length(ind.year)-z))
		res <- unlist(sapply(ind, function(z) apply(z, 1, do.budget)))
	}
	# output	
	return(res)
}