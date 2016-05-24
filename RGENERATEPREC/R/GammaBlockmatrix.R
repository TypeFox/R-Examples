# TODO: Add comment
# 
# Author: ecor
###############################################################################

NULL
#'
#' This return a \code{\link{blockmatrix}} object containing the gaussian cross-correlation matrices. 
#' 
#' @param data data frame or 'zoo' R object containing daily precipitation time series for several gauges (one gauge time series per column). See \code{\link{CCGamma}}.
#' @param lag numeric (expressed as number of days) used for the element [1,1] of the returned blockmatrix.  
#' @param p numeric order $p$ of the auto-regeression 
#' @param ... further argments of \code{\link{CCGamma}}
#' 
#' @details This a wrapper for \code{\link{CCGamma}} with the option \code{only.matrix=TRUE} and the function value is transformed into a \code{\link{blockmatrix}} object. 
#' 
#' @import blockmatrix Matrix
# @rdaname CCGamma
#' @seealso  \code{\link{CCGamma}},\code{\link[RMAWGEN]{continuity_ratio}},\code{\link{omega_inv}},\code{\link{omega}}
#' @examples 
#' 
#' data(trentino)
#' 
#' year_min <- 1961
#' year_max <- 1990
#' 
#' period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
#' station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
#' prec_mes <- PRECIPITATION[period,station]  
#' 
#' ## removing nonworking stations (e.g. time series with NA)
#' accepted <- array(TRUE,length(names(prec_mes)))
#' names(accepted) <- names(prec_mes)
#' for (it in names(prec_mes)) {
#' 		 accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
#' }
#'
#' prec_mes <- prec_mes[,accepted]
#' ## the dateset is reduced!!! 
#' prec_mes <- prec_mes[,1:2]
#' 
#' p <- 1 ## try p <- 2 !!!  
#' CCGamma <- CCGammaToBlockmatrix(data=prec_mes,lag=0,p=p,tolerance=0.001)
#' 
#' ## Not Run in the examples, uncomment to run the following line
#' # CCGamma_1 <- CCGammaToBlockmatrix(data=prec_mes,lag=1,p=p,tolerance=0.001)
#' 
#' 
#' ### Alternatively, recommended ..... 
#' ## Not Run in the examples, uncomment to run the following line
#' # CCGamma <- CCGammaToBlockmatrix(data=prec_mes,lag=0,p=p+1,tolerance=0.001)
#' 
#' # CCGamma0 <- CCGamma[1:p,1:p]
#' # CCGamma1 <- CCGamma[(1:p),(1:p)+1]
#' 
#' # CCGamma0_inv <- solve(CCGamma0)
#' 
#' 
#' ## Not Run in the examples, uncomment to run the following line
#' #a1 <- blockmatmult(CCGamma0,CCGamma0_inv)
#' # a2 <- blockmatmult(CCGamma1,CCGamma0_inv)
#' 
#' 
#' 
#' # CCGamma_1t <- t(CCGamma1)
#' #CCGamma_0t <- t(CCGamma0)
#' 
#' # A <- t(solve(CCGamma_0t,CCGamma_1t))
#' 
#' 
#' 
#' @export

CCGammaToBlockmatrix <- function(data,lag=0,p=3,...) {
	
	out <- NULL 
	

	matrix <- array(NA,c(p,p))
	for (r in 1:nrow(matrix)) {
		matrix[r,] <- lag[1]-r+(1:p)
	}
	
	
	
	lags <-   unique(as.vector(matrix)) # lag[1]+(-p+1):(p-1)
	out <- CCGamma(data,lag=lags,only.matrix=TRUE,...)
	
	if (!is.list(out)) {
		
		out <- list(out)
		names(out) <- lags[1]
	}
	
	names(out) <- paste("gamma",names(out),sep="_")
	matrix[,] <- paste("gamma",matrix,sep="_")
	

	
	out <- blockmatrix(value=matrix,list=out)
	
	return(out)
	
}