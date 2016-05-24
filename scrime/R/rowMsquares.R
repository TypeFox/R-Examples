`rowMsquares` <-
function(..., listTables=NULL, clScores=NULL, levScores=NULL, add.pval=TRUE){
	if(is.null(listTables))
		listTables <- list(...)
	if(!is.list(listTables))
		stop("listTables must be a list (of matrices).")
	n.class <- length(listTables)
	if(n.class<2)
		stop("At least two matrices need to be specified.")	
	cn <- checkTableList(listTables)
	n.lev <- length(cn)
	if("NA"%in%cn){
		warning("The column named NA is removed from all matrices.")
		ids <- which(cn=="NA")
		listTables <- lapply(listTables, function(x) x[,-ids])
	}
	if(is.null(levScores)){
		if(is.null(cn))
			stop("Either levScores must be specified or the matrices must have\n",
				"(the same) column names specifying numeric scores.")
		levScores <- as.numeric(cn)
		if(any(is.na(levScores)))
			stop("At least one of column names does not specify a numeric score.")
	}
	else{
		if(any(!is.numeric(levScores)))
			stop("levScores must be numeric.")
		if(length(levScores)!=n.lev)
			stop("The length of levScores must be equal to the number of columns.")
	}
	if(is.null(clScores))
		clScores <- 1:n.class
	else{
		if(length(clScores)!=n.class)
			stop("The length of clScores must be equal to the number of matrices.")
		if(any(!is.numeric(clScores)))
			stop("clScores must be numeric.")
	}	
	
	mat.n <- sapply(listTables, rowSums)
	n.obs <- rowSums(mat.n)
	ybar <- as.vector(mat.n%*%clScores) / n.obs
	X <- listTables[[1]]
	for(i in 2:n.class)
		X <- X + listTables[[i]]
	xbar <- as.vector(X%*%levScores) / n.obs
	mat.lev <- matrix(levScores, length(n.obs), n.lev, byrow=TRUE)
	mat.lev <- mat.lev - xbar
	listTables <- lapply(listTables, function(x) x * mat.lev)
	mat.class <- matrix(clScores, length(n.obs), n.class, byrow=TRUE)
	mat.class <- mat.class - ybar
	denom <- rowSums(mat.n * mat.class * mat.class)
	X <- listTables[[1]]
	num <- rowSums(X) * mat.class[,1]
	for(i in 2:n.class){
		X <- X + listTables[[i]]
		num <- num + rowSums(listTables[[i]]) * mat.class[,i]
	}
	denom <- denom * rowSums(X * mat.lev)
	stats <- num * num / denom
	stats <- stats * (n.obs-1)
	if(!add.pval)
		return(stats)
	rawp <- pchisq(stats, 1, lower.tail=FALSE)
	structure(list(stats=stats, rawp=rawp))
}

