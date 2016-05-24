# TODO: Add comment
# 
# Author: ecor
###############################################################################
NULL
#' Precipitation Occurence Model
#' 
#' This functions creates a stochastic Occurence Model for the variable  \code{x} (\code{PrecipitationOccurenceModel} S3 object) through a calibration from observed data.     
#' 
#' @param x variable utilized for the auto-regression of its occurence, e.g. daily precipitaton 
#' @param p auto-regression order 
#' @param exogen exogenous predictors
#' @param monthly.factor vector of factors indicating the month of the days
#' @param valmin minimum admitted value for daily precipitation amount 
#' @param id.name identification name of the station
#' @param ... further arguments
#' 
#' 
#' @return The function returns a \code{PrecipitationOccurenceModel-class} S3 object containing the following elements:
#' 
#' \code{predictor} data frame containg the endogenous and exogenous predictors of the logistic regression model;
#' 
#' \code{glm} the genaralized liner model using for the logistic regression;
#' 
#' \code{p} auto-regression order 
#' 
#' \code{valmin} minimum admitted value for daily precipitation amount 
#' 
#' @seealso \code{\link{glm}}
#' @export 
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
#' prec_occurence_mes <- prec_mes>=valmin
#' 
#' station <- names(prec_mes)[!(names(prec_mes) %in% c("day","month","year"))]

#' it <- station[2]
#' vect <- Tx_mes[,it]-Tn_mes[,it]
#' months <- factor(prec_mes$month)

#' model <- PrecipitationOccurenceModel(x=prec_mes[,it],exogen=vect,monthly.factor=months)
#' 
#'probs <- predict(model$glm,type="response")
#' 
#' 
#' 
#' 
#' 
#' plot(months[-1],probs)
#' 
#' newdata <- model$predictor[2000:2007,]
#' probs0 <- predict(model,newdata=newdata)
#' 
#' 
#' 

###

PrecipitationOccurenceModel <- function(x,exogen=NULL,p=1,monthly.factor=NULL,valmin=0.5,id.name=NULL,...) {
	
	if (is.null(valmin)) valmin <- NA
	if (!is.na(valmin)) variable <- x>=valmin
	
	
	out <- list()
	
	if (is.null(exogen) | (length(exogen)==0)) {
		
		df <- as.data.frame(array(variable,length(variable)))
		names(df) <- "variable"
		nrows <- nrow(df)
		
		
	##	out$exogen <- exogen
		
	} else  { 
	
		exogen <- as.data.frame(exogen)
		df <- exogen
		df$variable <- variable
		
		df <- df[,c("variable",names(exogen))]
		
		nrows <- nrow(df)
	##	out$exogen <- exogen
	}
	
	
	
	
	if (!is.null(monthly.factor)) df$month <- factor(monthly.factor)
	
	### lag analysis 
	out <- list()
	
	out$predictor <- df
	
	if (p>0) {
		ndf <- nrow(df)
		rows <- ((p+1):ndf)
		df <- df[rows,]
		for (l in 1:p) {
			
			label <- sprintf("x_l%02d",l)
			
			
			df[,label] <- variable[rows-l]
			out$predictor[,label] <- c(array(NA,p),variable[rows-l])
			
			
		}
		
	}
	
	###
	df <- data.frame(df)
	
	out$glm <- glm(df,family="binomial")
	
	out$p <- p 
	
	out$valmin <- valmin
##	 if (!is.null(monthly.factor)) out$month <- factor(monthly.factor)
	names_predictor <- names(out$predictor)
	out$predictor <- as.data.frame(out$predictor[,-1])
	names(out$predictor) <- names_predictor[-1]
	if (!is.null(id.name)) {
		
		out$id.name <-  id.name
		
	} else {
		
		out$id.name <-  NA
	}
	
	class(out) <- "PrecipitationOccurenceModel"
	return(out)
	
	
}

