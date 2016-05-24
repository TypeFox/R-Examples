`rowCATTs` <-
function(cases, controls, scores=NULL, add.pval=TRUE){
	if(!is.matrix(cases) | !is.matrix(controls))
		stop("cases and controls must be matrices.")
	if(any(dim(cases)!=dim(controls)))
		stop("cases and controls have not the same dimensions.")
	rn <- rownames(cases)
	if(any(rn!=rownames(controls)))
		stop("The row names differ between cases and controls.")
	cn <- colnames(cases)
	if(any(colnames(controls)!=cn))
		stop("The column names differ between cases and controls.")
	if("NA" %in% cn){
		ids <- which(cn=="NA")
		cases <- cases[,-ids]
		controls <- controls[,-ids]
		warning("The column named NA is removed from both cases and controls.")
	}
	if(is.null(scores)){
		if(is.null(cn))
			stop("Either scores must be specified or the matrices must have\n",
				"(the same) column names specifying numeric scores.")
		scores <- as.numeric(cn)
		if(any(is.na(scores)))
			stop("At least one of column names does not specify a numeric score.")
	}
	else{
		if(any(!is.numeric(scores)))
			stop("scores must be numeric.")
		if(length(scores)!=ncol(cases))
			stop("The number of scores must be equal to the number of columns.")
	}
	if(any(cases<0))
		stop("All values in cases must be non-negative integers.")
	if(any(controls<0))
		stop("All values in controls must be non-negative integers.")
	n.cases <- rowSums(cases)
	n.controls <- rowSums(controls)
	n <- n.cases + n.controls
	ybar <- n.cases / n
	xbar <- as.vector((cases + controls) %*% scores)
	xbar <- xbar / n
	mat.scores <- matrix(scores, length(n), length(scores), byrow=TRUE)
	mat.scores <- mat.scores - xbar
	cases <- cases * mat.scores
	controls <- controls * mat.scores
	num <- (-ybar) * rowSums(cases) + (1-ybar) * rowSums(controls)
	cases <- cases + controls
	denom <- rowSums(cases * mat.scores)
	denom <- denom * ybar * n.controls
	stats <- num * num / denom
	stats <- stats * n
	if(!add.pval)
		return(stats)
	rawp <- pchisq(stats, 1, lower.tail=FALSE)
	structure(list(stats=stats, rawp=rawp))
}

