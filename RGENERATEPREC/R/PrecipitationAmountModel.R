NULL
#' ....
#' 
#' @param x observed precipitation amount time series (data frame)
#' @param station string vector containing station identification codes
#' @param valmin maximum admitted value of precipitation depth 
#### @param method see \code{\link{cor}}
#' @param origin date of the day referred by he first row of \code{x}. 
#' @param sample character string. If it is \code{"monthly"} (Default), the corralaton matrix is calculeted per each month.  
#' @param ... further agruments for \code{\link[RMAWGEN]{normalizeGaussian_severalstations}}
#' 
#' 
#' @return The function returns AN S3 OBJECT ...... the correlation matrix of precipitation amount values (excluding the zeros). 
#' In case \code{sample=="monthly"} the runction return a \code{MonlthyList} S3 object.
#' 
#' @seealso \code{\link{predict.PrecipitationAmountModel}},\code{\link[RMAWGEN]{normalizeGaussian_severalstations}}
############# ,\code{\link{generate}},\code{\link{random.precipitation.values}},\code{\link{cor}}
#' 
#' @export
#' @importFrom RMAWGEN normalizeGaussian_severalstations
#' 
#'@examples 
#' library(RGENERATEPREC)
#' 
#' set.seed(1245)
#' 
#' data(trentino)
#' 
#' year_min <- 1961
#' year_max <- 1990
#' 
#' origin <- paste(year_min,1,1,sep="-")
#' end <- paste(year_max,12,31,sep="-")
#' 
#' period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
#' period_temp <- TEMPERATURE_MAX$year>=year_min & TEMPERATURE_MAX$year<=year_max
#' 
#' prec_mes <- PRECIPITATION[period,]
#' Tx_mes <- TEMPERATURE_MAX[period_temp,]
#' Tn_mes <- TEMPERATURE_MIN[period_temp,]
## removing nonworking stations (e.g. time series with NA)
#' accepted <- array(TRUE,length(names(prec_mes)))
#' names(accepted) <- names(prec_mes)
#' for (it in names(prec_mes)) {
#' 	acc <- TRUE
#' 	acc <- (length(which(!is.na(Tx_mes[,it])))==length(Tx_mes[,it]))
#' 	acc <- (length(which(!is.na(Tn_mes[,it])))==length(Tn_mes[,it])) & acc
#' 	accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it])) & acc
#' 	
#' }
#' 
#' valmin <- 1.0
###station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
#' prec_mes <- prec_mes[,accepted]
#' 
#' 
#' 
#' Tx_mes <- Tx_mes[,accepted]
#' Tn_mes <- Tn_mes[,accepted]
#' prec_occurence_mes <- prec_mes>=valmin
#' 
#' station <- names(prec_mes)[!(names(prec_mes) %in% c("day","month","year"))]
#' 
#' precamount <- PrecipitationAmountModel(prec_mes,station=station,origin=origin)
#' 
#' val <- predict(precamount)
#' 
#' prec_gen <- generate(precamount)  
#' 
#' 
#' 
#' 
#' month <- adddate(as.data.frame(residuals(precamount$T0090)),origin=origin)$month
#' #####plot(month,residuals(precamount$T0090))
#' plot(factor(month),residuals(precamount$T0090))
#' 
#' qqplot(prec_mes$T0083,prec_gen$T0083)
#' abline(0,1)
#' 
#' 
#' 
#' 





PrecipitationAmountModel <- function(x,valmin=1,station=names(x),sample="monthly",origin="1961-1-1",...) {
	
	if (!is.null(station)) {
		
		x <- x[,station]

		
	}
	amount <- x
	occurrence <- as.data.frame(x>=valmin)
	
	xm <- as.matrix(x)
	na.index     <- which(is.na(xm))
	dry.index  <- which(xm<valmin)
	
	

	
	xm[dry.index] <- NA
	
	x[,] <- xm
	
	gauss <- normalizeGaussian_severalstations(x=x,data=x,mean=0,sd=1,inverse=FALSE,sample=sample,origin_x=origin,origin_data=origin,...)
	
	
	
	
	#### FARE UN LAPPLY CON I MESI!!!!!
	
	
	
	
	
	
	if (sample=="monthly") {
		names <- names(occurrence)
		occurrence <- adddate(occurrence,origin=origin)
		month <- factor(occurrence$month)
		str(month)
		occurence <- occurrence[,names]
		occurence$month <- month
		
	}
	
##	names(station) <- station
##	out <- lapply(X=station,)
	out <- list()
	
	for ( it in station) {
				
				df <- occurence 
				df <- df[,names(df)!=it]
				df$gamount <- gauss[,it]
				
				
				names <- c("gamount",names(df)[names(df)!="gamount"])
				str(names)
				df <- df[,names]
				df <- df[!is.na(df$gamount),]
				
				out[[it]] <- lm(df)
				
				## set attribute with station
				attr(out[[it]],"station") <- it
				
				
				
	}			
 	#### 
	
	out$station <- station
	out$sample <- sample
	out$x <- amount
	out$origin <- origin
	out$valmin <- valmin
	class(out) <- "PrecipitationAmountModel"	
	return(out)
	

	
}


#NULL
#
##'
##' @param model 
##' 
##' 
##' @export 
##' 
##' 
#
#
#filling.precipitaation.amount <- function(model=NULL,occ,origin_occ,...) {
#	
#	
#	occ <- 
#	
#	
#	
#}













