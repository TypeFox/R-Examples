# ============================================================
# various association measures between sparse matrices
# ============================================================

# Pearson correlation matrix between columns of X, Y
# http://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r
#
# covmat uses E[(X-muX)'(Y-muY)] = E[X'Y] - muX'muY
# with sample correction n/(n-1) this leads to cov = ( X'Y - n*muX'muY ) / (n-1)
#
# the sd in the case Y!=NULL uses E[X-mu]^2 = E[X^2]-mu^2
# with sample correction n/(n-1) this leads to sd^2 = ( X^2 - n*mu^2 ) / (n-1)
#
# Note that results larger than 1e4 x 1e4 will become very slow, because the resulting matrix is not sparse anymore. 

corSparse <- function(X, Y = NULL, cov = FALSE) {

	X <- as(X,"dgCMatrix")
	n <- nrow(X)
	muX <- colMeans(X)
	
	if (!is.null(Y)) {
		stopifnot( nrow(X) == nrow(Y) )
		Y <- as(Y,"dgCMatrix")
		muY <- colMeans(Y)
		covmat <- ( as.matrix(crossprod(X,Y)) - n*tcrossprod(muX,muY) ) / (n-1)
		sdvecX <- sqrt( (colSums(X^2) - n*muX^2) / (n-1) )
		sdvecY <- sqrt( (colSums(Y^2) - n*muY^2) / (n-1) )
		cormat <- covmat/tcrossprod(sdvecX,sdvecY)
	} else {		
		covmat <- ( as.matrix(crossprod(X)) - n*tcrossprod(muX) ) / (n-1)
		sdvec <- sqrt(diag(covmat))
		cormat <- covmat/tcrossprod(sdvec)	
	}
	
	if (cov) {
		dimnames(covmat) <- NULL
		return(covmat)
	} else {
		dimnames(cormat) <- NULL
		return(cormat)
	}
}

# cosine similarity matrix between columns of X, Y
# results larger than 1e7 x 1e7 are just barely manageable on my laptop
# a random sparse 1e8 x 1e8 takes about 3 minutes, and 1.5 Gb Memory.
#
# different weightings of normalisations can be used	

	# allow for different weighting functions: can also be defined externally!
	# defined as a function of rowvalues (s) and number of columns (N)
	# no weighting (NULL) is taken as default
	idf <- function(s,N) { log(N/(1+s)) }
	# inverse square root
	isqrt <- function(s,N) { s^-0.5 }	
	# for use in later functions (e.g. sim.words)
	# weight = NULL leads to identical results, but is quicker
	none <- function(s,N) { s }
	
	# allow for different normalisation functions: can also be defined externally!
	# Euclidean 2-norm is taken as default
	norm2 <- function(x,s) { drop(crossprod(x^2,s)) ^ 0.5 }
	# Alternatively take 1-norm
	norm1 <- function(x,s) { drop(crossprod(abs(x),s)) }

cosSparse <- function(X, Y = NULL, norm = norm2 , weight = NULL) {

	X <- as(X,"dgCMatrix")
	if (!is.null(Y)) {
		stopifnot( nrow(X) == nrow(Y) )
		Y <- as(Y,"dgCMatrix")
	}

	if (!is.null(weight)) {	
		Nx <- ncol(X)
		Sx <- rowSums(abs(X))
		Wx <- Diagonal( x = match.fun(weight)(Sx,Nx) )
		X <- Wx %*% X
		if (!is.null(Y)) {
			Ny <- ncol(Y)
			Sy <- rowSums(abs(Y))
			Wy <- Diagonal( x = match.fun(weight)(Sy,Ny) )
			Y <- Wy %*% Y
		}
	}
			
	S <- rep(1,nrow(X))			
	N <- Diagonal( x = match.fun(norm)(X,S)^-1 )
	X <- X %*% N
	if (!is.null(Y)) {
		N <- Diagonal( x = match.fun(norm)(Y,S)^-1 )
		Y <- Y %*% N
		return(crossprod(X,Y))	
	} else {
		return(crossprod(X))
	}
}


cosMissing <- function(X, availX, Y = NULL, availY = NULL, norm = norm2 , weight = NULL) {

	X <- as(X,"dgCMatrix")
	if (!is.null(Y)) {
		stopifnot( nrow(X) == nrow(Y) )
		Y <- as(Y,"dgCMatrix")
	}

	if (!is.null(weight)) {	
		Nx <- ncol(X)
		Sx <- Nx * rowSums(abs(X)) / rowSums(availX)
		Wx <- Diagonal( x = match.fun(weight)(Sx,Nx) )
		X <- Wx %*% X
		if (!is.null(Y)) {
			Ny <- ncol(Y)
			Sy <- Ny * rowSums(abs(Y)) / rowSums(availY)
			Wy <- Diagonal( x = match.fun(weight)(Sy,Ny) )
			Y <- Wy %*% Y
		}
	}
	
	if (is.null(Y)) {
		N <- tcrossprod(match.fun(norm)(X,availX))
		R <- crossprod(X)
	} else {
		stopifnot(!is.null(availY))
		N <- tcrossprod(match.fun(norm)(X,availY),match.fun(norm)(Y,availX))
		R <- crossprod(X,Y)
	}
	N <- N*(as(R,"nMatrix")*1)
	R@x <- R@x/N@x
	return(R)
}

# rowGroup is a sparse matrix with same number of rows as X, or a vector of length nrow(X) specifying the grouping of the rows

cosRow <- function(X, rowGroup, Y = NULL, norm = norm2 , weight = NULL) {

	if (is.vector(rowGroup)) { 
		rowGroup <- t(ttMatrix(rowGroup)$M*1) 
	} else {
		rowGroup <- as(rowGroup,"dgCMatrix")
	}

	X <- as(X,"dgCMatrix")
	if (!is.null(Y)) {
		stopifnot( nrow(X) == nrow(Y) )
		Y <- as(Y,"dgCMatrix")
	}

	if (!is.null(weight)) {	
		N <- ncol(X)
		s <- rowSums(X)
		S <- N * s / ( drop(s %*% tcrossprod(rowGroup)) )
		W <- Diagonal( x = match.fun(weight)(S,N) )
		X <- W %*% X
		if (!is.null(Y)) {
			N <- ncol(Y)
			s <- rowSums(Y)
			S <- N * s / ( drop(s %*% tcrossprod(rowGroup)) )
			W <- Diagonal( x = match.fun(weight)(S,N) )
			Y <- W %*% Y
		}
	}
	
	if (is.null(Y)) {
		N <- tcrossprod(match.fun(norm)(X,rowGroup))
		R <- crossprod(X)
	} else {
		N <- tcrossprod(match.fun(norm)(X,rowGroup),match.fun(norm)(Y,rowGroup))
		R <- crossprod(X,Y)
	}
	N <- N*(as(R,"nMatrix")*1)
	R@x <- R@x/N@x
	return(R)
}

# colGroupX is either a sparse Matrix with the same number of columns as X, rows represent grouping. Or a grouping vector of length ncol(X), specifying group indices.

# weight does not really seem to make sense here: complete data each row should have an equal amount of entries! Any row-weighting would only measure the amount of missing data.

cosCol <- function(X, colGroupX, Y = NULL, colGroupY = NULL, norm = norm2 ) {
	
	X <- as(X,"dgCMatrix")
	if (is.vector(colGroupX)) {
		colGroupX <- ttMatrix(colGroupX)$M * 1 
	} else {
		colGroupX <- as(colGroupX,"dgCMatrix")
	}

	S <- rep(1,nrow(X))
	Freq <- tcrossprod(X,colGroupX)
	N <- Diagonal( x = drop(crossprod(colGroupX, match.fun(norm)(Freq,S)))^-1 )
	X <- X %*% N
	
	if (!is.null(Y)) {
		stopifnot( nrow(X) == nrow(Y) )
		stopifnot(!is.null(colGroupY))
		Y <- as(Y,"dgCMatrix")
		if (is.vector(colGroupY)) { 
			colGroupY <- ttMatrix(colGroupY)$M * 1 
		} else {
			colGroupX <- as(colGroupX,"dgCMatrix")
		}	

		Freq <- tcrossprod(Y,colGroupY)
		N <- Diagonal( x = drop(crossprod(colGroupY, match.fun(norm)(Freq,S)))^-1 )
		Y <- Y %*% N

		return(crossprod(X,Y))
		
	} else {
		return(crossprod(X))
	}	
}

# association matrix between columns of X, Y

	# allow for different functions to be specified: can be defined externally!
	# defaults to poisson
	poi <- function(o,e) { sign(o-e) * (o * log(o/e) - (o-e)) }
	# pointwise mutual information, aka "log-odds" in bioinformatics
	pmi <- function(o,e) { log(o/e) }
	# weighted pointwise mutual information, i.e the basis for MI
	wpmi <- function(o,e) { o * log(o/e)}
	# good old pearson residuals
	res <- function(o,e) { (o-e) / sqrt(e) }

# WATCH OUT! it might not always be the right decision to ignore the cases in which O=zero!!! e.g. residuals should not be zero then (but negative!!!).

assocSparse <- function(X, Y = NULL, method = res, N = nrow(X), sparse = TRUE) {

	X <- as(X,"ngCMatrix")*1
	Fx <- colSums(X)

	# observed coocurrences "O"
	if (is.null(Y)) {
		O <- crossprod(X)
		} else {
		stopifnot( nrow(X) == nrow(Y) )
		Y <- as(Y,"ngCMatrix")*1
		Fy <- colSums(Y)
		O <- crossprod(X,Y)
	}

	# The trick to keep things sparse here is to not compute anything when O is zero. In most practical situations this seems to be fine, but note that in cases in which the expectation is high, though observed is zero this might lead to large discrepancies.
	if (sparse) {	
		R <- as(O,"nMatrix")/N
		Fx <- Diagonal( x = Fx )
		if (is.null(Y)) {
			E <- Fx %*% R %*% Fx
			E <- as(E,"symmetricMatrix")
		} else {
			Fy <- Diagonal( x = Fy )
			E <- Fx %*% R %*% Fy
		}	
		R@x <- match.fun(method)(O@x,E@x)

	# this is the easy non-sparse method. Note that the result will be non-sparse, so this is not feasable for very large datasets
	} else {
		if (is.null(Y)) {
			E <- tcrossprod(Fx)/N
		} else {
			E <- tcrossprod(Fx,Fy)/N
		}
		R <- match.fun(method)(O,E)
	}
	return(R)
}

# still TO DO: assoc.missing

# The following is a version of assoc for cases in which the columns form groups. Typically found in case of a large collection of nominal variables (e.g. WALS) with an index matrix. If you want to establish the association between the variables, you'd need this version to get the expectation right.

# colGroupX is either a sparse Matrix with the same number of columns as X, rows represent grouping. Or a grouping vector of length ncol(X), specifying group indices.

assocCol <- function(X, colGroupX, Y = NULL, colGroupY = NULL, method = res, sparse = TRUE) {

	X <- as(X,"dgCMatrix")
	
	if (is.vector(colGroupX)) { 
		colGroupX <- ttMatrix(colGroupX)$M*1 
	} else {
		colGroupX <- as(colGroupX,"dgCMatrix")	
	}	
	Gx <- crossprod(colGroupX)
	
	# Fx: Frequencies of columns over all data, not just matched cases!
	# divided by 'N', occurrences of Fx, leading to probability of columns
	# will lead to slightly different values from traditional chisquare
	Fx <- colSums(abs(X))
	Fx <- Diagonal( x = Fx / drop(Gx %*% Fx) )

	if (is.null(Y)) {
		O <- crossprod(X)

		} else {
		stopifnot( nrow(X) == nrow(Y) )
		stopifnot(!is.null(colGroupY))

		Y <- as(Y,"dgCMatrix")
		
		if (is.vector(colGroupY)) { 
			colGroupY <- ttMatrix(colGroupY)$M*1
		} else {
			colGroupY <- as(colGroupY,"dgCMatrix")
		}		
		Gy <- crossprod(colGroupY)

		Fy <- colSums(abs(Y))
		Fy <- Diagonal( x = Fy / drop(Gy %*% Fy) )
		
		O <- crossprod(X, Y)
	}

	# in a sparse way: tricky to do, but possible. Ignore items where Observed is zero. Not that this is an approximation that might lead to somewhat unexpected results in specific situations.
	if (sparse) {	

		# Number of cases per block. To keep this sparse needs some tricks (unfold-refold)	
		if (is.null(Y)) {			
			N <- unfoldBlockMatrix(O, colGroupX)	
		} else {			
			N <- unfoldBlockMatrix(O, colGroupY)
		}
		G <- crossprod(kronecker( Diagonal(nrow(colGroupX)), colGroupX ))
		sums <- G %*% rowSums(N$U)
		S <- Diagonal( x = drop(sums) )
		N <- N$L %*% S %*% (as(N$U,"nMatrix")*1) # refold with 'weights' S

		# Expectation
		if (is.null(Y)) {
			E <- forceSymmetric(Fx %*% N %*% Fx)
		} else {
			E <- Fx %*% N %*% Fy
		}
		
		result <- O
		result@x <- match.fun(method)(O@x,E@x)	

	# in a non-sparse way: this is much more elegant! But the result is not sparse, so this is not feasable for very large datasets.
	} else {
		if (is.null(Y)) {
			E <- crossprod(X %*% Gx %*% Fx)
		} else {
			E <- crossprod(X %*% Gx %*% Fx, X %*% Gy %*% Fy)
		}
		result <- match.fun(method)(O,E)
	}
	return(result)
}

# same as above here, but for groups of rows. E.g. with WALS when looking at the similarity between languages, given the expectation of groups (features) as 1/nr.of.values 

# rowGroup is a sparse matrix with same number of rows as X, or a vector of length nrow(X) specifying the grouping of the rows

assocRow <- function(X, rowGroup, Y = NULL, method = res) {

	X <- as(X,"dgCMatrix")
	
	if (is.vector(rowGroup)) {
		rowGroup <- t(ttMatrix(rowGroup)$M*1)
	} else {
		rowGroup <- as(rowGroup,"dgCMatrix")
	}
	W <- Diagonal( x = rowSums(tcrossprod(rowGroup))^-0.5 )

	if (is.null(Y)) {
		O <- crossprod(X) 	
		E <- crossprod(t(rowGroup) %*% W %*% X)
	} else {
		stopifnot( nrow(X) == nrow(Y) )
		Y <- as(Y,"dgCMatrix")
		O <- crossprod(X,Y)
		E <- crossprod(t(rowGroup) %*% W %*% X, t(rowGroup) %*% W %*% Y)
	}

	E <- E * (as(O,"nMatrix")*1) # this is necessary for cases in which O is zero
	O@x <- match.fun(method)(O@x,E@x)
	return(O)
}
