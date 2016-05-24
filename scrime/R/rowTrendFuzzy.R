rowTrendFuzzy <-
function(score, probs, y, mat.fuzzy=NULL, alternative=c("two.sided", "less", "greater"),
		check=TRUE){
	if(is.null(mat.fuzzy)){
		if(missing(score) | missing(probs))
			stop("score and probs need to be specified, if mat.fuzzy is not specified.")
		mat.fuzzy <- getMatFuzzy(score, probs, check=check)
		rm(probs)
	}
	else{
		if(!missing(score) | !missing(probs))
			stop("If score and probs are specified, mat.fuzzy is not allowed to be specified.")
		if(!is.matrix(mat.fuzzy))
			stop("mat.fuzzy must be a matrix.")
		if(any(mat.fuzzy < 0))
			stop("All values in mat.fuzzy must be non-negative.")
	}
	if(any(is.na(mat.fuzzy)))
		stop("No missing values allowed in mat.fuzzy (or in probs).")
	if(any(is.na(y)))
		stop("No missing values allowed in y.")
	if(any(!y %in% c(0,1)))
		stop("The values in y must be either 0 or 1.")
	if(length(y) != ncol(mat.fuzzy))
		stop("The length of y differs from the number of columns of mat.fuzzy (or in probs).")
	piHat <- mean(y)
	num <- as.vector(mat.fuzzy %*% (y-piHat))
	mat.fuzzy <- mat.fuzzy - rowMeans(mat.fuzzy)
	denom <- rowSums(mat.fuzzy * mat.fuzzy) * piHat * (1-piHat)
	z <- num / sqrt(denom)
	alt <- match.arg(alternative)
	rawp <- switch(alt, less=pnorm(z), greater=pnorm(z, lower.tail=FALSE), 
		two.sided=2*pnorm(-abs(z)))
	names(num) <- names(z)
	structure(list(stat=z, rawp=rawp, theta=num/denom, varTheta=1/denom))
}

