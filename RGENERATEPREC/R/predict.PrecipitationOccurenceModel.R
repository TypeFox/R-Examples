
NULL

#' @title Prediction of a  \code{PrecipitationOccurenceModel} model object
#' 
#' @description It is a wrapper of \code{\link{predict.glm}} method for the a \code{PrecipitationOccurenceModel} model object S3 class. 
#' 
#' @param object model returned by \code{\link{PrecipitationOccurenceModel}}
#' @param newdata predictor or exogenous variables
#' @param type see \code{\link{predict.glm}}. Default is \code{"response"}. See \code{\link{predict.glm}}.
#' @param previous logical vector containing previously occurred states.
#' @param endogenous String vector containing the name of the endogenous variables. 
#' It is used if the endogenous variables are more than one, otherwise is set \code{NULL}(Default).
#' @param ... further arguments 
#' 
#' @seealso \code{\link{predict.glm}},\code{\link{PrecipitationOccurenceModel}}
#' @export 
#' @method predict PrecipitationOccurenceModel
#' @S3method predict PrecipitationOccurenceModel
#' @aliases predict predict.PrecipitationOccurenceModel 
#' 
#' @rdname predict
##### @importFrom predict stats
#' @importFrom RMAWGEN normalizeGaussian_severalstations
#' 
#' @seealso \code{\link{predict.glm}},\code{\link{predict.glm}},\code{\link{PrecipitationOccurenceModel}},\code{\link{PrecipitationAmountModel}}
#' 
#' @examples
#' 
#' library(RGENERATEPREC)
#' 
#' data(trentino)
#' 
#' year_min <- 1961
#' year_max <- 1990
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
#' origin <- paste(year_min,1,1,sep="-")
#' 
#' 
#' prec_occurence_mes <- prec_mes>=valmin
#' 
#' station <- names(prec_mes)[!(names(prec_mes) %in% c("day","month","year"))]

#' it <- station[2]
#' vect <- Tx_mes[,it]-Tn_mes[,it]
#' months <- factor(prec_mes$month)

#' model <- PrecipitationOccurenceModel(x=prec_mes[,it],exogen=vect,monthly.factor=months)
#' 
#' probs <- predict(model)
#' 
#' nday <- 3.0
#' vect_new <- array(1.0,nday)
#' months_new <- array(1,nday)
#' row_test <- 2000:2007
#' newdata <- model$predictor[row_test,]
#' probs2 <- predict(model,newdata=newdata)
#' 
#' probs[row_test]==probs2
#' ###
#' 
#' 
#' 
#' 
#'
#' prec_occurence_mes <- prec_mes>=valmin
#' 
#' station <- names(prec_mes)[!(names(prec_mes) %in% c("day","month","year"))]
#' 
#' station <- station[1:4] ## reduced the dataset!!! 
#' Tx_mes <- Tx_mes[,station]
#' Tn_mes <- Tn_mes[,station]
#' 
#' prec_mes <- prec_mes[,station]
#' exogen <- Tx_mes-Tn_mes
#' months <- factor(prec_mes$month)
#' 
#' 
#' ### Not Run 
#' ### Please uncomment the following lines to run them 
#' 
#' 
#' #model_multisite <- PrecipitationOccurenceMultiSiteModel(x=prec_mes,
#' #exogen=exogen,origin=origin,multisite_type="wilks")
#' #
#' #
#' #model_multisite_logit <- PrecipitationOccurenceMultiSiteModel(x=prec_mes,
#' #exogen=exogen,origin=origin,multisite_type="logit")
#' #
#' #
#' #probs_multimodel  <- predict(model_multisite_logit)
#' 
#' 


predict.PrecipitationOccurenceModel <- function(object,newdata=NULL,type="response",previous=NULL,endogenous=NULL,...) {
	
	if (object$p<1) previous <- NULL
	
	
	if (!is.null(previous)) {
		
		previous <- as.vector(previous)
		
		previous[!is.logical(previous)] <- previous[!is.logical(previous)]>=object$valmin
		
		if(is.null(newdata)) newdata <- object$predictor
		
		newdata <- as.data.frame(newdata)
		
		if (is.null(endogenous)) {
			for (l in 1:object$p) {
			
			label <- sprintf("x_l%02d",l)
			
			newdata[,label] <- previous[l]
			
			}
		
			
		
			labels <- sprintf("x_l%02d",1:object$p)
		
			names(newdata)[!(names(newdata) %in% labels)] <- names(object$predictor)[!(names(object$predictor) %in% labels)]
		} else { 
			labels <- NULL
			for (l in 1:object$p) {
				
				######rows_l <- ((l+1):ndf)
				
				label <- sprintf("_endog_l%02d",l)
				names_label <- paste(endogenous,label,sep="")
				newdata[,names_label] <- previous[,l]
				labels <- c(labels,names_label)
			}	
		
			names(newdata)[!(names(newdata) %in% labels)] <- names(object$predictor)[!(names(object$predictor) %in% labels)]
		}	
		
	} else {
		
		if(is.null(newdata)) newdata <- object$predictor
		
		newdata <- as.data.frame(newdata)
		
		names(newdata) <- names(object$predictor)
		
	}
	
	
	newdata <- newdata[,names(newdata) %in% names(object$predictor)]
#####	newdata <- newdata[,names(object$predictor)]
	
	out <- NULL
	out <- predict(object$glm,newdata=newdata,type=type,...)
	if (is.null(newdata) & (object$p>0)) {
		
		out <- c(array(NA,object$p),out)
		
	}
	return(out)
}


NULL

#'
#' 
#' 
#' @export 
#' @method predict PrecipitationOccurenceMultiSiteModel
#' @S3method predict PrecipitationOccurenceMultiSiteModel
#' @aliases predict predict.PrecipitationOccurenceMultiSiteModel
#' @rdname predict
#'
predict.PrecipitationOccurenceMultiSiteModel <- function(object,...) { ## exogen=NULL commented
	   	
	out <- NULL
	if (object$type=="wilks") {
		endogenous <- NULL
	} else if (object$type=="logit") {
	
	 	endogenous <- object$station
	
	}
	
	out <- lapply(X=object[object$station],FUN=predict,endogenous=endogenous,...)
	names(out) <- object$station
	
	out <- as.data.frame(do.call(cbind,out))
	
	names(out) <- object$station
	
	
	return(out)
	
	
	
	
	
}
