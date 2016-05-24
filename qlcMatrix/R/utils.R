# ============================================================
# Special "reduced" KhatriRao version in which empty rows are removed.
# ============================================================

# This is only important for the new rownames, as now only rownames are produced for the non-empty rows. That is more efficient than making all rownames, as in the offical version in the Matrix package. It is *extremely* tricky to get the names right: watch out with the order of the indices!

rKhatriRao <- function(X, Y
                       , rownamesX = rownames(X)
                       , rownamesY = rownames(Y)
                       , simplify = FALSE
                       , binder = ":"
                       , FUN = "*"
                       ) {

	# sparse KhatriRao
	M <- KhatriRao(X, Y, FUN = FUN)

	# remove empty rows
	selection <- rowSums(M, sparseResult = TRUE) > 0
	M <- M[selection@i,]
	
	# make names for the non-empty rows
	nonzero <- Matrix(	selection,
						nrow = nrow(Y),
						ncol = nrow(X),
						sparse = TRUE
						)
	nonzero <- as(nonzero,"TsparseMatrix")

	rownamesM <- paste(	rownamesY[nonzero@i + 1],
						rownamesX[nonzero@j + 1],	
						sep = binder
						)
	
	if (simplify) {
		rownames(M) <- rownamesM
		colnames(M) <- colnames(X)
		return(M)
	} else {
		return(	list(	M = M,
						rownames = rownamesM
						))
	}
}

# ====================================================
# Construct random sparse matrices, useful for testing
# ====================================================

# code from Martin Maechler

rSparseMatrix <- function(nrow, ncol, nnz, 
                          rand.x = function(n) round(rnorm(nnz), 2), ...)
{
    stopifnot((nnz <- as.integer(nnz)) >= 0,
              nrow >= 0, ncol >= 0, nnz <= nrow * ncol)
    if (is.null(rand.x)) {
    	sparseMatrix(i = sample(nrow, nnz, replace = TRUE),
                     j = sample(ncol, nnz, replace = TRUE),
                     dims=c(nrow,ncol))
    } else {
    sparseMatrix(i = sample(nrow, nnz, replace = TRUE),
                 j = sample(ncol, nnz, replace = TRUE),
                 x = rand.x(nnz), dims = c(nrow, ncol), ...)
    }
}

# =====================================================
# Unfold blockmatrix
# first by column groups, optionally also by row groups
# =====================================================

unfoldBlockMatrix <- function(X, colGroups, rowGroups = NULL) {

	if (is.vector(colGroups)) {
		colGroups <- ttMatrix(colGroups)$M
	} else {
		colGroups <- as(colGroups, "dgCMatrix")
	}

	U <- KhatriRao(colGroups,X)
	L <- as(kronecker( t(rep(1,nrow(colGroups)))
	                   , Diagonal(nrow(X)) )
	        ,"CsparseMatrix")
		
	if (is.null(rowGroups)) {
		return( list( U=U, L=L ) )
	} else {
		
		if (is.vector(rowGroups)) {
			rowGroups <- ttMatrix(rowGroups)$M
		} else {
			rowGroups <- as(rowGroups, "dgCMatrix")
		}

		R <- as(kronecker( rep(1,nrow(rowGroups))
		                   , Diagonal(ncol(U)))
		        ,"CsparseMatrix")
		
		rowGroups <- kronecker( t(rep(1,nrow(colGroups)))
		                        , rowGroups )
		U <- t( KhatriRao(rowGroups , t(U)))

		return( list( U=U, L=L, R=R ))
	}
}

# ======================
# Maximum per row/column
# ======================

# returns sparse vector with maximum values.
# Optionally returns a sparse matrix 
# with the position of these maxima in the original matrix
# becomes very slow when number of entries in the table is larger than 1e5.

rowMax <- function(X, which = FALSE, ignore.zero = TRUE) {

# old approach, much slower
# new approach much faster using "rollup" in package "slam" !!!
#	m <- aggregate(x~i, data = summary(X), FUN = max)
#	maximum <- sparseVector(x = m$x, i = m$i, length = nrow(X))

	Y <- as.simple_triplet_matrix(drop0(X))
	maximum <- as(slam::rollup(Y, 2, FUN = max), "sparseVector")
	
	if(!ignore.zero) {
		maximum@x[maximum@x<0] <- 0
	}
	
	if (which) {
		d <- Diagonal(x = as(maximum,"vector"))
		W <- as(X,"nMatrix") * 1
		Xmax <- d %*% W
		W@x <- (X@x == Xmax@x) * 1
		W <- as(drop0(W), "nMatrix")
		return(list(max = maximum, which = W))
		} else {
			return(maximum)
	}
}

rowMin <- function(X, which = FALSE, ignore.zero = TRUE) {

# old approach, much slower
#	m <- aggregate(x~i, data = summary(X), FUN = min)
#	minimum <- sparseVector(x = m$x, i = m$i, length = nrow(X))

	Y <- as.simple_triplet_matrix(drop0(X))
	minimum <- as(slam::rollup(Y, 2, FUN = min), "sparseVector")
	
	if(!ignore.zero) {
		minimum@x[minimum@x > 0] <- 0
	}
	
	if (which) {
		d <- Diagonal(x = as(minimum, "vector"))
		W <- as(X, "nMatrix") * 1
		Xmin <- d %*% W
		W@x <- (X@x == Xmin@x) * 1
		W <- as(drop0(W), "nMatrix")
		return(list(min = minimum, which = W))
		} else {
			return(minimum)
	}
}

colMax <- function(X, which = FALSE, ignore.zero = TRUE) {
	tmp <- rowMax(t(X), which = which, ignore.zero = ignore.zero)
	if (which) {
		return(list(max = tmp$max, which = t(tmp$which)))
	} else {
		return(tmp)
	}
}

colMin <- function(X, which = FALSE, ignore.zero = TRUE) {
	tmp <- rowMin(t(X), which = which, ignore.zero = ignore.zero)
	if (which) {
		return(list(min = tmp$min, which = t(tmp$which)))
	} else {
		return(tmp)
	}
}
