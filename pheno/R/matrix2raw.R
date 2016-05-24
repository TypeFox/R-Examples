# Converts a matrix M into a dataframe D with three columns (x, factor 1, factor 2)
# where rows of M are ranks of factor 1 levels and columns of M are
# ranks of factor 2 levels, missing values are assumed to be 0 or NA.
# Returned data frame as ordered first by factor 2 and then factor 1.
# input: matrix M with two optional ordered numeric vectors for level names.
matrix2raw <- function(M,l1,l2) {
	if(!is.matrix(M)) stop("matrix2raw: first argument must be a matrix. Exiting ...")
	if(missing(l1)) { l1 <- as.vector(c(1:dim(M)[1]),"numeric") }
	else { l1 <- as.vector(l1,"numeric") }
	if(!is.vector(l1,"numeric") || length(l1) != dim(M)[1]) {
		stop("matric2raw: check argument l1, not numeric or wrong length. Exiting ...")
	}
	if(missing(l2)) { l2 <- as.vector(c(1:dim(M)[2]),"numeric") }
	else { l2 <- as.vector(l2,"numeric") }
	if(!is.vector(l2,"numeric") || length(l2) != dim(M)[2])
		stop("matric2raw: check argument l2, not numeric or wrong length. Exiting ...")

	n <- 0
	y <- f1 <- f2 <- vector("numeric",dim(M)[1]*dim(M)[2])
	for(j in 1:dim(M)[2]) {			# order first by factor 2 
		for(i in 1:dim(M)[1]) {		# then by factor 1
			if(is.na(M[i,j])) { M[i,j] <- 0 }
			if(M[i,j] != 0) {
				n <- n + 1 
				y[n] <- M[i,j]; f1[n] <- l1[i]; f2[n] <- l2[j]
			}
		}
	}
	y <- y[1:n]; f1 <- f1[1:n]; f2 <- f2[1:n]
	D <- data.frame(y,f1,f2)
	return(D)
}
