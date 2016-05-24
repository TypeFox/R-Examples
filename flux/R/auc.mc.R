auc.mc <- function(x, y, method = "leave out", lo = 2, it = 100, ...){
	METHODS <- c("leave out", "bootstrap", "sorted bootstrap", "constrained bootstrap", "jackknife", "jack-validate")
	method <- pmatch(method, METHODS)
	if (is.na(method)){
		stop("invalid method")
		}
	
	if(method == 1){
		ind <- lapply(seq(it), function(z) sort(sample(seq(length(x)), length(x)-lo)))
		res <- sapply(ind, function(z) auc(x[z], y[z], ...))
	}
	else if (method == 2){
		ind <- lapply(seq(it), function(z) sample(seq(length(x)), replace=TRUE))
		res <- sapply(ind, function(z) auc(x, y[z], ...))
	}
	else if (method == 3){
		ind <- lapply(seq(it), function(z) sort(sample(seq(length(x)), replace=TRUE)))
		res <- sapply(ind, function(z) auc(x, y[z], ...))
	}
	else if (method == 4){
		ind <- lapply(seq(it), function(z) unique(sort(sample(seq(length(x)), replace=TRUE))))
		res <- sapply(ind, function(z) auc(x[z], y[z], ...))
	}	
	else if (method == 5){
		ind <- combs(seq(length(x)), length(x)-1)
		res <- apply(ind, 1, function(z) auc(x[z], y[z], ...))
	}	
	else if (method == 6){
		ind <- lapply(c(1:lo), function(z) combs(seq(length(x)), length(x)-z))
		res <- unlist(sapply(ind, function(z) apply(z, 1, function(w) auc(x[w], y[w], ...))))
	}		
	
	return(res)
}