NULL
#'
#' This function calculates  ...
#' 
#' @param data data frame or 'zoo' R object containing daily precipitation time series for several gauges (one gauge time series per column). See \code{\link{CCGamma}}.
#' @param CCGamma0 correlation block-matrix with lag of 0 days. Object returned by \code{\link{CCGammaToBlockmatrix}}. If omitted,default is \code{NULL}, it is internally calculated.
#' @param CCGamma1 correlation block-matrix with lag of 1 days. Object returned by \code{\link{CCGammaToBlockmatrix}}. If omitted,default is \code{NULL}, it is internally calculated.
#' @param p numeric order $p$ of the auto-regeression, see \code{\link{CCGammaToBlockmatrix}}
#' @param sample character string indicated if the coefficients must be estimated differently for subset of the year, e.g. monthly. Admitted values are \code{NULL} (Default), \code{"all"} or \code{"monthly"}. 
#' @param origin character string (yyyy-dd-mm) indicated the date of the first row of \code{"data"}. It is used if \code{data} and \code{sample} are not \code{NULL}.
#' @param ... other arguments of \code{\link{CCGammaToBlockmatrix}}
#' 
#' 
#' @author Emanuele Cordano
#' 
#' @return A S3 object of class \code{"YuleWalkerCoefficientBlockmatrices"} (or \code{"YuleWalkerCoefficientBlockmatricesPerEachMonth"} in case \code{sample="monthly"}) which is a list containing the block matrices \code{A},\code{Sigma_u} of the Yule-Walker Equation and the object \code{CCGammaInfo} containing probabilities of no precipitation occurence and returned by function \code{\link{CCGamma}} applied with \code{lag=0}. In case \code{sample="monthly"}) 
#' functioion return a \code{"YuleWalkerCoefficientBlockmatricesPerEachMonth"}, i. e. a list of \code{"YuleWalkerCoefficientBlockmatrices"} for each month. 
#' 
#' @note This function uses Yule-Walker equations for VAR to estimate the coefficient block-matrices blockmatrix \code{A} and \code{Sigma_u}. The input of this function are the correletion block-matrices \code{CCGamma0} and \code{CCGamma1}. 
#' If they are missing (and then \code{NULL}) , they are also calculated from the original dataset (argument \code{data}). In this last case, the coefficients can be estiomated differently for each monthly setting \code{sample} equal to \code{"monthly"}.
#' 
#' 
#' @seealso \code{\link{CCGammaToBlockmatrix}},\code{\link{generatePrecipitationAmount}}
#' 
#' @export
#' @examples 
#' 
#' library(RMRAINGEN)
#' 
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
#' # ## Not Run in the examples, uncomment to run the following line
#' # coeff <- CoeffYWeq(data=prec_mes,p=1,tolerance=0.001)
#' 
#' #
#' # 
#' # Alternatively the coefficients of Vector Auto-Regressive Model 
#' # can be separately calculated for each month   
#' 
#' # ## Not Run in the examples, uncomment to run the following line
#' #origin <- paste(year_min,1,1,sep="-")
#' #
#' #
#' 
#' #coeff_monthly <- CoeffYWeq(data=prec_mes,p=1,tolerance=0.001,sample="monthly",origin=origin)
#' 
#'

# 
#dmvnorm(x=c(0,0))
#dmvnorm(x=c(0,0), mean=c(1,1))
#
#sigma <- matrix(c(4,2,2,3), ncol=2)
#x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma)
#colMeans(x)
#var(x)
#
#x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma, method="chol")
#colMeans(x)
#var(x)
#
#plot(x)
#

CoeffYWeq <- function(data=NULL,CCGamma0=NULL,CCGamma1=NULL,p=1,sample=NULL,origin="1961-01-01",...) {
	
	
##
#   sampling 
#
#
#
##	

	format="%Y-%m-%d"

	if (is.null(sample)) sample <- "all"

	if (sample=="monthly" & !is.null(data)) {
		
		if (is.na(origin)) {
			
			origin <- as.character(index(date)[1],format=format)
			
			
		}
		
		data <- as.data.frame(data)
		names <- names(data)
		print(origin)
		data <- adddate(data, origin = origin)
		
	#	origin0 <- as.character(as.POSIXct(origin,tz="A")-3600*24*p,format=format)
	#	data0 <- adddate(data, origin = origin0)
		months <- 1:12
	   	out <- lapply(X=months,FUN=function(m,data,p,names,...) {
				   
				   mprev <- m-1
				   mnext <- m+1
				   mprev[mprev==0] <- 12
				   mnext[mnext==13] <- 1
					
				   cond <- data$month==m
				   cond <- ((data$month==mprev) & (data$day>=25)) | cond
	   			   cond <- ((data$month==mnext) & (data$day<=8)) | cond	   
				   data <- data[cond,]
				   data <- data[,names]
				   print(m)
				#   print(data)
				 #  str(data)
				   out <- CoeffYWeq(data=data,p=p,sample="all",origin=NA,...)	
		##         
		
				  return(out)
		
		
		
		},data=data,p=p,names=names,...)
		
		names(out) <- months
		class(out) <- "YuleWalkerCoefficientBlockmatricesPerEachMonth"
 		return(out)
		
		## END EXECUTION IN MONTHLY OPTION
		
	}
	
	
	
	
	
	
	
	out <- list(A=NULL,Sigma_u=NULL)
	
	 
	
	if ((is.null(CCGamma0)) | (is.null(CCGamma1))) {
		
		
		
		CCGamma <- CCGammaToBlockmatrix(data,lag=0,p=p+1,...) 
		
		CCGamma0 <- as.blockmatrix(CCGamma[1:p,1:p],nrow=p,ncol=p)
		CCGamma1 <- as.blockmatrix(CCGamma[(1:p),(1:p)+1],nrow=p,ncol=p)
###		CCGamma2 <- as.blockmatrix(CCGamma[(1:p),(1:p)+2],nrow=p,ncol=p)
		
		
	}

####	CCGamma_2t <- t(CCGamma2)
	CCGamma_1t <- t(CCGamma1)
    CCGamma_0t <- t(CCGamma0)
  
	
	
 


	out$A <- t(solve(CCGamma_0t,CCGamma_1t,symm.precond=TRUE))

	A <- as.matrix(out$A)
	eigen <- eigen(A)
	amax <- max(abs(eigen$values))
	THRS <- 0.99999
	if (amax>=0.9999) {
		warning("Eigenvalues of A are greater than 1 and then put less than 1 for VAR stability ( A is modified)!!!")
		ev <- eigen$values
		ev[abs(eigen$values)>amax] <- amax
		A <- eigen$vectors %*% diag(ev) %*% solve(eigen$vectors) 
		out$A <- as.blockmatrix(A,nrow=p,ncol=p)
		
	}


##	if (max(eigens)>=1) stop("Eigrnvalues greater than 1 in mod!!")
	
	
	out$Sigma_u <- CCGamma0-blockmatmult(out$A,CCGamma_1t)
	
	
	###### rainfall occurence 
	
	
	out$CCGammaInfo <- NULL
	
	if (!is.null(data)) {
		
		out$CCGammaInfo <- CCGamma(data,lag=0,p=NA,only.matrix=FALSE,...)
		
		
	}
	
	
	
	class(out) <- "YuleWalkerCoefficientBlockmatrices"
	return(out)

	
	## here is the blockmatrix
	
	
	
	
	
	
}
