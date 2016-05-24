# Converts a numeric data frame D with three columns (x, factor 1, factor 2)
# in a matrix M where rows are ranks of levels of factor 1 and columns are
# ranks of levels of factor 2, missing values are set to NA.
raw2matrix <- function(D) {
	if(!is.data.frame(D) || length(D) != 3) 
		stop("raw2matrix: argument must be data frame with 3 fields. Exiting ...")

	n <- length(rownames(D))
	f1 <- factor(D[[2]])
	nrows <- nlevels(f1)
	f2 <- factor(D[[3]])
	ncols <- nlevels(f2)
	M <- matrix(nrow=nrows,ncol=ncols)
	for( k in 1:n) {
		i <- f1[[k]]; j <- f2[[k]]; M[i,j] <- D[k,1]
	}
	return(M)
}
