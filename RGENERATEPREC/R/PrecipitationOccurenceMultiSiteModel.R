
#
NULL
#' Precipitation Occurence Multi-Site Model
#' 
#' This functions creates a stochastic Occurence Multi-Site Model for the variable  \code{x} (\code{PrecipitationOccurenceMultiSiteModel} S3 object) through a calibration from observed data.     
#' 
#' @param x data frame (each column is a site) of variable utilized for the auto-regression of its occurrence, e.g. daily precipitaton 
#' @param exogen exogenous predictors
#' @param station character string vectors containing the codes of the station used for model calibration
#' @param valmin minimum admitted value for daily precipitation amount 
#' @param multisite_type string indicating the utilized approach for spatial multi-site dependence description. Default is \code{"wilks"}.
#' @param tolerance_wilks see \code{tolerance} used by \code{\link{omega_inv}} through \code{\link{CCGamma}}
#' @param origin character string (yyyy-dd-mm) indicating the date of the first row of \code{"x"}.
#' @param p auto-regression order 
#' @param ... further arguments
#' 
#' 
#' @return The function returns a \code{PrecipitationOccurenceModel-class} S3 object containing the following elements:
#' 
#' ... \code{\link{PrecipitationOccurenceModel}} S3 class objects for each analyzed site. The name is the site (or station) code
#' 
#' \code{ccgama} \code{CCGammaObjectListPerEachMonth} object, i.e. matices of Gaussian Inter-Site  Correlation returned by \code{\link{CCGamma}};
#' 
#' \code{type} string indicating the utilized approach for spatial multi-site dependence description, only \code{"wilks"} type is implemented;
#' 
#' \code{station} character string vectors containing the codes of the station used in \code{PrecipitationMultiSiteOccurenceModel}.
#' 
#'  @export
#' 
#' @seealso \code{\link{PrecipitationOccurenceModel}},\code{\link{CCGamma}} 
#' 
#' @import stringr
#' @examples
#' 
#' library(RGENERATEPREC)
#' 
#' data(trentino)
#' 
#' year_min <- 1961
#' year_max <- 1990
#' origin <- paste(year_min,1,1,sep="-")
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

#' station <- station[1:2] # to save example elapsed time!!
#' exogen <- Tx_mes-Tn_mes
#' months <- factor(prec_mes$month)
#' #' ### Not Run!! 
#' # The following lines are commented to save example elapsed time!!
#' #model_multisite <- PrecipitationOccurenceMultiSiteModel(x=prec_mes,exogen=exogen,
#' #origin=origin,multisite_type="wilks")
#' 
#' ### Not Run!! 
#' #  The following lines are commented to save example elapsed time!!
#' #model_multisite_logit <- PrecipitationOccurenceMultiSiteModel(x=prec_mes,exogen=exogen,
#' #origin=origin,multisite_type="logit")
#' ### 
PrecipitationOccurenceMultiSiteModel <- function(x,exogen=NULL,station=names(x),origin=origin,valmin=0.5,multisite_type="wilks",tolerance_wilks=0.001,p=2,...) {
	
	x <- adddate(x,origin=origin)
	monthly.factor <- factor(x$month)
	station <- station
	station <- station[!(station %in% c("day","month","year"))]
	
	x <- x[,station]
	
	if (multisite_type=="wilks") {
		
	
		if (is.null(exogen) | (length(exogen)==0)) {
			
			exogen <- lapply(X=station,FUN=function(x){ NULL })
			names(exogen) <- station
			
			
		}
	
		out <- list()
	
		
		for (it in station) {
			
			if (is.data.frame(exogen)) {
				
				cols <- str_detect(names(exogen),it)
				exogen_loc <- exogen[,cols]
				
				
			} else if (is.list(exogen)) {
				
				exogen_loc <- exogen[[it]]
			} else {
				
				exogen_loc <- exogen
			}
			
			out[[it]] <- PrecipitationOccurenceModel(x=x[,it],exogen=exogen_loc,monthly.factor=monthly.factor,valmin=valmin,p=p,id.name=it,...)
			
		}
		
#		CCGamma(data, lag = 0, p0_v1 = NULL, p = NA, valmin = 0.5,
# 				nearPD = (lag >= 0), interval = c(-1, 1),
# 				tolerance = .Machine$double.eps, only.matrix = FALSE,
# 				return.value = NULL, null.gcorrelation = 1e-05, sample = NULL,
# 				origin = "1961-1-1", ...)
#
#
#		
#
		
		out$ccgamma <- CCGamma(x,lag=0,tolerance=tolerance_wilks,sample="monthly",only.matrix=TRUE,origin=origin,valmin=valmin)
		
		names(out$ccgamma) <- sprintf("month%02d",1:length(out$ccgamma))
	} else if (multisite_type=="logit"){
		
		if (is.null(valmin)) valmin <- NA
		if (!is.na(valmin))  { 
			variable <- x>=valmin
			x <- as.data.frame(variable)
		}
		
		if (!is.null(exogen)) { as.data.frame(array(0*NA,c(nrow(x),0))) }
		if (!is.data.frame(exogen)) {
			
			stop("Option: multisty_type == logit, exogen must be a data frame or NULL!!!")
		} else if (nrow(x)!=nrow(exogen)) {
			
			stop("Option: multisty_type == logit, exogen and x must be have the same numbers of rows!!!")
			
		}
		
		df <- exogen 
	###	str(x)
	###	str(df)
	###	str(p)
		### TO GO ON MONDAY ......
		if (p>0) {
			ndf <- nrow(df)
			rows <- ((p+1):ndf)
	
			for (l in 1:p) {
				rows_l <- ((l+1):ndf)
				
				label <- sprintf("_endog_l%02d",l)
				names_label <- paste(station,label,sep="")
		
				df[,names_label]  <- NA
				df[rows_l,names_label] <- x[rows_l-l,station]
	
	
	
	#######			out$predictor[,label] <- c(array(NA,p),variable[rows-l])
				
				
			}
		}
		
		
		out <- list()
		for (it in station) { 
			
			out[[it]] <- PrecipitationOccurenceModel(x=x[,it],exogen=df,monthly.factor=monthly.factor,valmin=valmin,p=0,...)
			out[[it]]$p <- p
		}
		
		#str(out)
		#stop()
		
#		} if (is.null(exogen) {
#			
#		} else if (!)
#		
#		out
		
		
		
		
		
	} else {
		
		out <- list()
	}
	
	
	out$K <- ncol(x)
	out$type <- multisite_type
	out$station <- station 
	out$p <- p 
	class(out) <- "PrecipitationOccurenceMultiSiteModel"
	return(out)
	
	
	
}

