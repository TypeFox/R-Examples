# TODO: Add comment
# 
# Author: ecor
###############################################################################

NULL
#'
#' This function calculates the parameters of a parametric probability distribution by fitting daily precipitation data for each station and each month. It is a wrapper function for \code{\link{fitdistr}} 
#' 
#' @param data dataset
#' @param dname name of the pobability distribution to be fitted, e.g. \code{"exp"}. 
#' @param start initialization configuration of probability distribution parameters for \code{\link{fitdistr}}.
#' @param sample character string indicated if the parameters  must be estimated differently for subset of the year, e.g. monthly. Admitted values are \code{NULL} (Default), \code{"all"} or \code{"monthly"}. 
#' @param origin character string containing the date releted to the first row of \code{data}. Default is \code{"1961-1-1"}.
#' @param valmin threshold precipitation value [mm] for wet/dry day indicator. 
#' If precipitation is lower than \code{valmin}, day is considered dry. Default is 0.5 mm. See \code{\link{continuity_ratio}},\code{\link{CCGamma}}.
#' @param ... further arguments for \code{\link{fitdistr}}
#' 
#' @return a list containig the fitting parameters: S3  \code{fitdistrForEachStation-class} or \code{fitdistrForEachStationForEachMonth-class} (\code{sample=="monthly"}) object
#' 
#' 
#' 
#' @seealso \code{\link{fitdistr}},\code{\link{adddate}}
#' 
#' @export
#' 
#' 
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
#' prec_mes <- prec_mes[,1:3]
#' 
#' fit  <- fitdistrForPrecipitation(data=prec_mes,dname="exp",start=NULL,sample=NULL) 
#' 
#' origin <- paste(year_min,1,1,sep="-")
#' 
#' 
#' fit_monthly  <- fitdistrForPrecipitation(data=prec_mes,dname="exp",
#' 					start=NULL,sample="monthly",origin=origin)
#' fit_monthly_gamma  <- fitdistrForPrecipitation(data=prec_mes,dname="gamma",
#' 						 start=NULL,sample="monthly",origin=origin)
#'
#' 
#' 
#' 
#' 
fitdistrForPrecipitation <- function(data,dname="exp",start=NULL,sample="all",origin="1961-1-1",valmin=0.5,...) {
	
	out <- NULL
	dfun <- get(paste("d",dname,sep=""))
	pfun <- get(paste("p",dname,sep=""))


 	if(is.null(sample)) sample <- "all"
	
	if (sample=="all") {
		
		if (is.null(start)) {
			start <- formals(dfun)
			start <- start[!(names(start) %in% c("x","log"))]
			
		}
		
		if (dname=="gamma") { 
			
			if (!is.numeric(start$rate)) {
				
				start <- start[!(names(start) %in% c("rate"))]
				
			} else if (!is.numeric(start$shape)) {
				
				start <- start[!(names(start) %in% c("scale"))]
				
			} else { 
				
				
				remove <- -which(names(start) %in% c("scale","rate"))[2]
				
				start <- start[remove]
			}
		}
		start <- lapply(X=start,FUN=function(x){ 
		if (!is.numeric(x)) {
							
			x <- 1
		}
		return(x)
		})
		
	    data <- as.data.frame(data)
		out <- list()
		for (c in 1:ncol(data)) {
			
			x <- data[,c]
			x <- x[x>valmin]
			
			out[[c]] <- fitdistr(x=x,dfun,start=start,...)
			out[[c]]$name <- dname 
			
		}
		names(out) <- names(data)
		class(out) <- "fitdistrForEachStation"
		return(out)
		
	} else if (sample=="monthly") {
		
		data <- adddate(data,origin=origin)
		ignore.date <- !(names(data) %in% c("year","month","day"))
		out <- list()
		months <- sort(unique(data$month))
		for (m in months) {
			
			val <- data[which(data$month==m),ignore.date]
			
			out[[m]] <- fitdistrForPrecipitation(val,dname=dname,start=start,sample="all",valmin=valmin,...) 
			
			
		}
		names(out) <- months
		class(out) <- "fitdistrForEachStationForEachMonth"
		return(out)
		
		
		
	}
	
	
	
	
	return(out)
}

#
#distributionPar <- lapply(X=names(distributions),FUN=function(x) {
#			dfun <- get(paste("d",x,sep=""))
#			out <- formals(dfun)
#			out <- out[!(names(out) %in% c("x","log"))]
#			for (f in 1:length(out)) {
#				
#				if (class(out[[f]])!="numeric") out[[f]] <- 1
#			}
#			
#			out$densfun <- dfun
#			out$namefun <- x
#			
#			return(out)
#		})
#names(distributionPar) <- names(distributions)
#
#
### Specific correction on density function arguments 
#distributionPar$gamma <- distributionPar$gamma[names(distributionPar$gamma)!="scale"]
#
#
#
#distributionParForEachStation <- list()
#distributionGofForEachStation <- list()
#for (it in stations) {
#	data <- prec_mes[,it]
#	data <- data[data>0]
#	## fitting parameters
#	distributionParForEachStation[[it]] <-  lapply(X=distributionPar,FUN=function(x,data,signif=0.05){
#				
#				out <- list()
#				start <- x[!(names(x) %in% c("densfun","namefun"))]
#				pnamefun <- paste("p",x$namefun,sep="")
#				print(pnamefun)
#				out$fit <- fitdistr(data,x$densfun,start=start) ##,method="mle") ###,method="CG") ###,fix.arg=args.fix)
#				
#				out$gof <- gof.test.mod(x=data,y=pnamefun,par=out$fit$estimate,what="ks.test")
#				
#				
#				return(out)
#				
#			},data=data)				
#	
#	
#	
#}
