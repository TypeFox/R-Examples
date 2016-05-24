"sos" <- function(x, method="mean", foc=NULL, normal.sp=TRUE, normal.pl=TRUE) {
	x <- x[,colSums(x)>0]
	subtra <- switch(method,
		mean = apply(x, 2, mean, na.rm=TRUE),
		foc = as.numeric(x[foc,]))
	res <- sum((t(x)-subtra)^2)
	if(normal.sp){ res <- res/ncol(x) }
	if(normal.pl){ res <- res/nrow(x) }
	names(res) <- "sososm"
	res
	}