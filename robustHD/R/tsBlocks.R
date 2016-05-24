# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Construct predictor blocks for time series models
#' 
#' Construct blocks of original and lagged values for autoregressive time 
#' series models with exogenous inputs.  The typical use case is to supply the 
#' output as \code{newdata} argument to the 
#' \code{\link[=predict.tslars]{predict}} method of robust groupwise least 
#' angle regression models.
#' 
#' @param x  a numeric matrix or data frame containing the exogenous predictor 
#' series.
#' @param y  a numeric vector containing the response series.
#' @param p  an integer giving the number of lags to include (defaults to 2).
#' @param subset  a logical or integer vector defining a subset of observations 
#' from which to construct the matrix of predictor blocks.
#' @param intercept  a logical indicating whether a column of ones should be 
#' added to the matrix of predictor blocks to account for the intercept.
#' 
#' @return A matrix containing blocks of original and lagged values of the 
#' time series \code{y} and \code{x}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{predict.tslars}}, \code{\link{tslars}}, 
#' \code{\link{predict.tslarsP}}, \code{\link{tslarsP}}
#' 
#' @keywords ts
#' 
#' @export

tsBlocks <- function(x, y, p = 2, subset = NULL, intercept = TRUE) {
    # initializations
    n <- length(y)
    x <- as.matrix(x)
    if(nrow(x) != n) stop(sprintf("'x' must have %d rows", n))
    if(!is.numeric(p) || length(p) == 0) p <- 2
    p <- as.integer(p[1])
    if(p < 1) {
        p <- 1
        warning("lag length too small, using lag length 1")
    }
	# take subset (if supplied)
	if(!is.null(subset)) {
		x <- x[subset, , drop=FALSE]
		y <- y[subset]
    	n <- length(y)
	}
    # add response to predictor matrix
    x <- cbind(y, addColnames(x))
    cn <- colnames(x)
    # check if we only have one row in the output and need to correct for that
    if(n-p < 0) stop("not enough observations")
    correct <- n-p == 0
    if(correct) rn <- rownames(x)[n]
    # construct blocks and combine them columnwise
    x <- lapply(seq_len(ncol(x)), function(i) block(x[, i], p))
    if(correct) {
        # make sure blocks are row vectors
        x <- lapply(x, t)
    }
    x <- do.call(cbind, x)
    if(correct) rownames(x) <- rn
    # set column names
    colnames(x) <- blockNames(cn, p)
    # add intercept
    if(isTRUE(intercept)) x <- addIntercept(x)
    # return matrix
    x
}


# construct a block of original and lagged values from time series
# x ... time series
# p ... lag length
block <- function(x, p) {
    n <- length(x)
    sapply(seq_len(p), function(ip) x[(p-ip+1):(n-ip+1)])
}

# default names for time series blocks
# cn ... vector of column names
# p .... lag length
blockNames <- function(cn, p) {
    if(p == 1) cn
    else c(rbind(cn, sapply(cn, paste, seq_len(p-1), sep=".")))
}
