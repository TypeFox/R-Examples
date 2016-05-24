#' Transformation of Factors to Individual Levels
#' 
#' The function \code{cv.trans.psa} takes a covariate data frame and replaces each 
#' categorical covariate of \code{n >=3} levels with n new binary covariate columns, 
#' one for each level. Transforms covariate dataframe for use with the function 
#' \code{cv.bal.psa}.
#' 
#' NOTE: This function originated in the \code{PSAgraphics} package. It has been
#' adapted here for the \code{multilevelPSA} package.
#' 
#' @param covariates A dataframe of covariates, presumably some factors.
#' @param fcol An optional vector containing the factor columns in the covariate 
#'        data frame. In \code{NULL} (default) routine to identfy factors internally.
#' 
#' @author James E. Helmreich James.Helmreich@@Marist.edu
#' @author Robert M. Pruzek RMPruzek@@yahoo.com
#' @author KuangNan Xiong harryxkn@@yahoo.com
#' @author Jason Bryer jason@@bryer.org
cv.trans.psa <- function (covariates, fcol = NULL) {
	X <- covariates
	n.rows <- dim(X)[1]
	col2 <- fcol
	if (is.null(fcol)) {
		xclass <- sapply(X, is.numeric)
		for (i in 1:dim(X)[2]) {
			if (!xclass[i]) {
				col2[i] <- i
			}
			else {
				col2[i] <- 0
			}
		}
	}
	col2 <- col2[col2[] != 0]
	if (sum(col2) == 0) {
		stop("No Categorical/Factor Columns Identified")
	}
	fac.num <- length(col2)
	d <- matrix(0, nrow = dim(X)[1], ncol = fac.num)
	fac.size <- NULL
	fac.codes <- NULL
	for (i in 1:fac.num) {
		fac.i <- X[, col2[i]]
		fac.levels <- matrix(sort(unique(fac.i)))
		num.levels <- length(fac.levels)
		colnames(fac.levels) <- dimnames(X)[[2]][col2[i]]
		rownames(fac.levels) <- paste(dimnames(X)[[2]][col2[i]], 
									  1:num.levels)
		fac.codes <- rbind(fac.codes, fac.levels)
		colnames(fac.codes) <- c("Levels")
		if (num.levels == 2) {
			fac.size[i] = 1
		} else {
			fac.size[i] = num.levels
		}
		for (j in 1:n.rows) {
			for (k in 1:num.levels) {
				if(is.na(fac.levels[k]) | is.na(fac.i[j])) {
					d[j, i] <- NA
				} else if(fac.levels[k] == fac.i[j]) {
					d[j, i] = k - 1
				}
			}
		}
	}
	dd = matrix(0, nrow = n.rows, ncol = sum(fac.size))
	dlabel = matrix(0, nrow = 1, ncol = sum(fac.size))
	ko = 0
	for (p in 1:fac.num) {
		h1 = ko + 1
		ko = ko + fac.size[p]
		h2 = ko
		if (fac.size[p] > 2) {
			fac.p <- X[, col2[p]]
			mmdp <- matrix(0, n.rows, length(unique(d[, p])))
			dd[, h1:h2] <- ifelse(d[, p] == (col(mmdp) - 1), 
								  1, 0)
			dlabel[, h1:h2] = paste(dimnames(X)[[2]][col2[p]], 
									sort(unique(fac.p)), sep = "_")
		}
		if (fac.size[p] == 1) {
			fac.p <- X[, col2[p]]
			dd[, h2] = d[, p]
			dlabel[, h2] = paste(dimnames(X)[[2]][col2[p]], sort(unique(fac.p))[2], 
								 sep = "_")
		}
		X[, col2[p]] = d[, p]
	}
	colnames(dd) = dlabel[1, ]
	X2 = cbind(X, dd)
	X = X2[, -col2]
	out <- list(X)
	names(out) <- c("covariates.transformed")
	return(out)
}
