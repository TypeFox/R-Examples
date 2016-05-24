# ------------------------------------------------------------------------------------------------------------
# >>>>>        Multi-dimensional MA normalization for plate effect        <<<<<
# ------------------------------------------------------------------------------------------------------------
# Normalize data to minimize the difference between sample plates (batch effects)
# - the primary method is Multi-MA, but other fitting function, f in manuscript (e.g. loess) is available, too.
# ------------------------------------------------------------------------------------------------------------

normn_MA <- function(mD, expGroup, represent_FUN= function(x) mean(x, na.rm= T), fitting_FUN= NULL, isLog= TRUE) {
	
	stopifnot(!missing(mD))
	stopifnot(!missing(expGroup))
	stopifnot(nrow(mD) == length(expGroup))
	if(!inherits(expGroup, "factor")) expGroup <- factor(expGroup)
	
	if(isLog) mD <- log(mD)

	# matrix of representative values of every experimental group 
	represent_FUN <- match.fun(represent_FUN)
	X <- sapply(levels(expGroup), function(x) apply(mD[expGroup == x, ], 2, represent_FUN))
	
	## NA full column -> excluded in normalization
	naCol <- apply(X, 2, function(ij) all(is.na(ij)))
	if(any(naCol)) {
		warning(paste(names(naCol)[naCol], "contains only NAs. The warning above is likely due to that."))
		X <- X[, !naCol]
	}
	
	nP <- ncol(X)	# the number of the experimental groups (e.g. plates)
	stopifnot(nP > 1)	# This function is designed to normalize multiple dimensional data
	
	# A = x1 + x2 + x3 + x4 + ... / nP = (the vector that pass through origin and line of identity)
	
	# >> Orth : an Orthogonal matrix for projection onto the subspace perpendicular to A <<
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Orth <- t( svd(rep(1, nP), nu = nP)$u )		# incl. A = u[, 1]
	# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if(Orth[1,1] < 0) Orth <- -1*Orth			# let the first element positive

	# >>> M and A values <<<
	# ~~~~~~~~~~~~~~~~~~~~~~~~
	MA <- t( Orth %*% t(X) )			# the values of X on M and A coordinates
	M <- as.matrix(MA[ ,2:nP])			# matrix
	A <- MA[ , 1]						# vector
	# ~~~~~~~~~~~~~~~~~~~~~~~~

	## M(normalized) = M(raw) - rbind(0, "Md")
	##                                    **
	Md <- sapply(1:(nP-1), function(j) {
		m_j <- M[, j]
		# >>> Find fitted value <<<
		# ****************************************************
		if(is.null(fitting_FUN)) (m_j)
		else {
			fitting_FUN <- match.fun(fitting_FUN)
			fitting_FUN(m_j= m_j, A= A)
		}
		# ****************************************************
	})
	Md <- rbind(0, t(Md))
	
	## M = (Orth) %*% X
	## X(normalized) = X - solve(Orth) %*% Md    (solve(Orth) = t(Orth))
	##                     ******************
	Xn <- t(t(Orth) %*% Md)

	## NA full column -> NA column
	if(any(naCol)) {
		for(ii in c(which(naCol) - 1)) {
			Xn <- if(ii == ncol(Xn)) {
				cbind(Xn, matrix(nrow= nrow(Xn), ncol= 1))
			} else {
				cbind(Xn[, 1:ii], matrix(nrow= nrow(Xn), ncol= 1), Xn[, (ii+1):ncol(Xn)])
			}
		}
	}
	
	colnames(Xn) <- colnames(X)

	# normalize individual values
	x_normn <- mD - t(Xn)[expGroup, ]
	if(isLog) x_normn <- exp(x_normn)
	
	return(invisible( x_normn ))
}
