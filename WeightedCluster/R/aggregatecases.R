wcAggregateCases <- function(x, weights=NULL, ...){
	UseMethod("wcAggregateCases")
}

wcAggregateCasesInternal <- function(x, weights=NULL){
	x <- as.data.frame(x)
	lx <- nrow(x)
	if(is.null(weights)){
		weights <- rep(1, lx)
	}
	for(i in 1:ncol(x)){
		x[, i] <- factor(x[ ,i])
		levels(x[, i]) <- as.character(1:nlevels(x[, i])) 
	}
	ids <- apply(x, 1, paste, collapse="@@@WC_SEP@@")
	FuncEnv <- environment()
	mcorr <- rep(NA, lx)
	myfunction <- function(x){
		FuncEnv$mcorr[x] <- x[1]
		return(c(x[1], sum(weights[x])))
	}
	xx <- aggregate(1:lx, list(id=ids), myfunction)$x
	mcorr2 <- match(mcorr, xx[ ,1])
	ret <- list(aggIndex=xx[, 1], aggWeights=xx[, 2], disaggIndex=mcorr2, disaggWeights=weights)
	class(ret) <- c("wcAggregateCases", class(ret))
	return(ret)
}

wcAggregateCases.default <- function(x, weights=NULL, ...){
	xx <- wcAggregateCasesInternal(x, weights=weights)
	return(xx)
}


wcAggregateCases.data.frame <- function(x, weights=NULL, ...){
	xx <- wcAggregateCasesInternal(x, weights=weights)
	return(xx)
}
wcAggregateCases.matrix <- function(x, weights=NULL, ...){
	xx <- wcAggregateCasesInternal(x, weights=weights)
	return(xx)
}

wcAggregateCases.stslist <- function(x, weights=NULL, weighted=TRUE, ...){
	if(is.null(weights) && weighted) {
		weights <- attr(x, "weights")
	}
	xx <- wcAggregateCasesInternal(x, weights=weights)
	return(xx)
}

print.wcAggregateCases <- function(x, ...){
	cat("Number of disaggregated cases: ", length(x$disaggWeights), "\n")
	cat("Number of aggregated cases: ", length(x$aggWeights), "\n")
	cat("Average aggregated cases: ", format(length(x$disaggWeights)/length(x$aggWeights)), "\n")
	cat("Average (weighted) aggregation: ", mean(x$aggWeights), "\n")
}
