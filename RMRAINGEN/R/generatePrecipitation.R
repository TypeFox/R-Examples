# TODO: Add comment
# 
# Author: ecor
###############################################################################


#'
#' Generetion of daily precipitation amount taking into accont a previous generation of precipitation occurence with \code{\link{generate}}
#' 
#' @title Generation of Precipitation Amount 
#' 
#' @param x object of precipitation quantiles returned by \code{\link{generate.YuleWalkerCoefficientBlockmatrices}} or \code{\link{generate.YuleWalkerCoefficientBlockmatricesPerEachMonth}} with option \code{precipitation.indicator==TRUE} 
#' @param origin character string containing the date releted to the first row of \code{x}. Default is \code{"1961-1-1"}.
#' @param par \code{fitdistrForEachStationForEachMonth-class} or \code{fitdistrForEachStation-class} object returned by \code{\link{fitdistrForPrecipitation}}
#' @param ... further arguments
#' @export
#' 
#' 
#' 
#' 
#' @seealso \code{\link{qqfit}},\code{\link{fitdistrForPrecipitation}},\code{\link{generate}}
#' 
#' 
#' @examples 
#' 
#' 
#' ## # Not Run in the examples, uncomment to run the following lines
#' 
#' #library(RMRAINGEN)
#' #
#' #set.seed(1233)
#' #data(trentino)
#' #
#' #year_min <- 1987
#' #year_max <- 1990
#' #
#' #period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
#' #station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
#' #prec_mes <- PRECIPITATION[period,station]  
#' #
#' ### removing nonworking stations (e.g. time series with NA)
#' #accepted <- array(TRUE,length(names(prec_mes)))
#' #names(accepted) <- names(prec_mes)
#' #for (it in names(prec_mes)) {
#' #		 accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
#' #}
#'
#' #prec_mes <- prec_mes[,accepted]
#' ### the dateset is reduced!!! 
#' #prec_mes <- prec_mes[,1:2]
#' #
#' #fit  <- fitdistrForPrecipitation(data=prec_mes,dname="exp",start=NULL,sample=NULL) 
#' #
#' #origin <- paste(year_min,1,1,sep="-")
#' #
#' ### Fitting of Probability Distribution of Precipitation Amount
#' #fit_monthly   <- fitdistrForPrecipitation(data=prec_mes,dname="gamma",
#' #                 start=NULL,sample="monthly",origin=origin)
#' #
#' ### Estimate coefficients for Precipitation Occurence Modeling  
#' ### (using generate.YuleWalkerCoefficientBlockmatrices)
#' #coeff_monthly <- CoeffYWeq(data=prec_mes,p=1,tolerance=0.001,sample="monthly",origin=origin)
#' #
#' #generation_monthly <- generate(coeff_monthly,year_min=year_min,year_max=year_max,
#' #					  names=names(prec_mes),precipitation.indicator=TRUE)
#' # 
#' #prec_gen <- generatePrecipitationAmount(x=generation_monthly,origin=origin,par=fit_monthly)
#' #
#' ### Estimate coefficients for Precipitation Occurence Modeling  (using generate.CCGammaObject)
#' #CCGamma_monthly <- CCGamma(data=prec_mes,lag=0,tolerance=0.001,only.matrix=FALSE,
#' #					sample="monthly",origin=origin)
#' #
#' #generation_monthly_2 <- generate(CCGamma_monthly,year_min=year_min,year_max=year_max,
#' #						names=names(prec_mes),precipitation.indicator=TRUE)
#' # 
#' #prec_gen_2 <- generatePrecipitationAmount(x=generation_monthly_2,origin=origin,par=fit_monthly)
#' #
#' ### Check Q-Q plots between observations and generations
#' #idst <- names(prec_mes)[2]
#' #month <- 6 ## June
#' #momth <- c(12,1,2) # Winter
#' #
#' #prec_mes <- adddate(prec_mes,origin=origin)
#' #prec_gen <- adddate(prec_gen,origin=origin)
#' #prec_gen_2 <- adddate(prec_gen_2,origin=origin)
#' #
#' #
#' #qqplot(prec_mes[prec_mes$month %in% month,idst],prec_gen[prec_gen$month %in% month,idst])
#' #abline(0,1)
#' #qqplot(prec_mes[prec_mes$month %in% month,idst],prec_gen_2[prec_gen_2$month %in% month,idst]) 
#' #abline(0,1)
#' #
#' ### Not Run in the examples, uncomment to run the following line
#' ### CCGamma_monthly_gen <- CCGamma(data=prec_gen,lag=0,tolerance=0.001,
#' ###                                only.matrix=FALSE,sample="monthly",origin=origin)
#' ### CCGamma_monthly_gen_2 <- CCGamma(data=prec_gen_2,lag=0,tolerance=0.001,
#' ###                                only.matrix=FALSE,sample="monthly",origin=origin)
#' #
#' #









generatePrecipitationAmount <- function(x,origin="1961-1-1",par=NULL,...){
	
	out <- NULL
	x <- as.data.frame(x)
	stations <- names(x)
	out <- x
	
	
	
	
	if (class(par)=="fitdistrForEachStation") {
		
		for (it in stations) {
			
		####	FUN <- get(paste("q",par[[it]]$name,sep=""))
		####	out[,it] <- FUN()
			str(x)
			str(out)
			print(it)
			out[,it] <- qqfit(as.vector(x[,it]),FUN=par[[it]]$name,par=par[[it]]$estimate,use.x=TRUE)
			
		}
		
	} else if (class(par)=="fitdistrForEachStationForEachMonth") {
		x <- adddate(x,origin=origin)
		months <- unique(x$month)
		
		for ( m in months) {
			
			out[x$month==m,] <- generatePrecipitationAmount(x=x[x$month==m,stations],origin=NULL,par=par[[m]])
			
			
		}
		
	}
	
	return(out)
}

