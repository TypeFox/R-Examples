NULL
#'
#'  Calculates the continuity ratio of a set of precipitation measured or generated data in several sites as defined by Wilks, 1998 (see reference link)
#' 
#' @param data containing daily precipitation time series for several gauges (one gauge time series per column)
#' @param lag numeric lag (expressed as number of days) used for computation for "cross" continuity ratio and joint probability of prercipitation (no)occurence. 
#' @param valmin threshold precipitation value [mm] for wet/dry day indicator. 
#' If precipitation is lower than \code{valmin}, day is considered dry. Default is 0.5 mm.
#' 
#' @author Emanuele Cordano, Emanuele Eccel 
#' 
#' @export 
#' 
#' @note If \code{lag==0} the function returns the continuity ratio and joint probability as described by Wilks, 1998. Otherwise the precipitation values for each couple of rain gauges are taken with \code{lag}-day lag. 
#' 
#' @references  see the following URL references:  \url{http://onlinelibrary.wiley.com/doi/10.1002/joc.2305/abstract} 
#'   and \url{http://www.sciencedirect.com/science/article/pii/S0022169498001863}
#' @return A list containing the following matrices: 
#' 
#' \code{continuity_ratio} : \code{lag}-day lagged  continuity ratio , 
#'
#' \code{occurence}  : joint probability of \code{lag}-day lagged precipitation occurence  
#' 
#' \code{nooccurence} : joint probability of \code{lag}-day lagged no precipitation occurence.
#' 
#' \code{noccurence_occurence} : joint probability of \code{lag}-day lagged no precipitation and precipitation occurence respectively.
#' 
#' \code{occurence_nooccurence} : joint probability of \code{lag}-day lagged precipitation and no precipitation occurence respectively.
#' 
#' \code{probability_continuity_ratio}: \code{lag}-day lagged ratio about precipitation probability contitioned to no precipitation/preciitation occurence in the other site
#' 
#' 
#' 
#' @examples 
#' 
#' data(trentino)
#' 
#' year_min <- 1961
#' year_max <- 1990
#' origin <- paste(year_min,1,1,sep="-")
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
#' continuity_ratio <-continuity_ratio(data=prec_mes,lag=0,valmin=0.5)
#' 
#' 
#' 





continuity_ratio <- function(data,lag=0,valmin=0.5) {
	
	ncols <- ncol(data)
	nrows <- nrow(data)
	out <- new.env()
	out$continuity_ratio <- array(NA,c(ncols,ncols))
	out$probability_continuity_ratio <- array(NA,c(ncols,ncols))
	out$occurence <- array(NA,c(ncols,ncols))
	out$nooccurence <- array(NA,c(ncols,ncols))
	out$nooccurence_occurence <- array(NA,c(ncols,ncols))
	out$occurence_nooccurence <- array(NA,c(ncols,ncols))
	lagp <- abs(lag)
	for (i in 1:ncols) {
		for (j in 1:ncols) {
			
			d1 <- data[(lagp+1):nrows,i]
			d2 <- data[1:(nrows-lagp),j]
			
			e1 <- mean(d1[which((d2>=valmin) & (d1>=valmin))],na.rm=TRUE)## e1 <- mean(d1[d1>valmin & d2>valmin],na.rm=TRUE)
			e2 <- mean(d1[which((d2<valmin) & (d1>=valmin))],na.rm=TRUE) ##e2 <- mean(d1[d1>valmin & d2<=valmin],na.rm=TRUE) ## something is wrong here!!!!
			l1 <- length(d1[which(d2>=valmin)])  ##length(d1[d1>valmin & d2>valmin])
			l2 <- length(d1[which(d2<valmin)]) ##length(d1[d1>valmin & d2<=valmin])
			
			el1 <- mean(as.vector(d1[which(d2>=valmin)])>=valmin,na.rm=TRUE)
			el2 <- mean(as.vector(d1[which(d2<valmin)])>=valmin,na.rm=TRUE)
			nrowsa <- length(d1[!is.na(d1) & !is.na(d2)])
			
			out$continuity_ratio[i,j] <- e2/e1
			out$probability_continuity_ratio[i,j] <- el2/el1
			out$occurence[i,j] <- length(d1[d1>=valmin & d2>=valmin & !is.na(d1) & !is.na(d2)])/nrowsa
			out$nooccurence[i,j] <- length(d1[d1<valmin & d2<valmin & !is.na(d1) & !is.na(d2)])/nrowsa
			out$nooccurence_occurence[i,j] <- length(d1[d1<valmin & d2>=valmin & !is.na(d1) & !is.na(d2)])/nrowsa
			out$occurence_nooccurence[i,j] <- length(d1[d1>=valmin & d2<valmin & !is.na(d1) & !is.na(d2)])/nrowsa
		}
		
	}
	if (lag<0) {
		
		out$continuity_ratio <- t(out$continuity_ratio)
		out$probability_continuity_ratio <- t(out$probability_continuity_ratio)
		out$occurence <- t(out$occurence)
		out$nooccurence <- t(out$nooccurence)
		out$noccurence_occurence <- t(out$noccurence_occurence)
		out$occurence_noccurence <- t(out$occurence_noccurence)
		
	} 
	return(as.list(out))
}


