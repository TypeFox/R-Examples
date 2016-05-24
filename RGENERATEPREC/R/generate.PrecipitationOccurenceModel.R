
NULL 

#' Stochastic Generation of a  \code{PrecipitationOccurenceModel} or \code{PrecipitationOccurenceMultiSiteModel} model object
#' 
#' It is an implentation of \code{\link{generate}} method 
#' 
#' @param x model returned by \code{\link{PrecipitationOccurenceModel}} or \code{\link{PrecipitationOccurenceMultiSiteModel}}
#' @param newdata predictor or exogenous variables. See \code{\link{predict.PrecipitationOccurenceModel}}
#' @param exogen  predictor or exogenous variables
#' @param monthly.factor vector of factors indicating the month of the days
#' @param random vector of random or calculated numbers ranging between 0 and 1 
#' @param origin,end  character strings (yyyy-dd-mm) indicating the start and/or end date of the daily weather generation.
#' @param n number of generations. See \code{\link{generate}}. Here it is ignored and the number of generations is given by \code{origin},\code{end} or \code{monthly.factor}.
#' @param previous logical vector containing previously occurred states
#' @param ... further arguments 
#' 
#' @seealso \code{\link{generate}},\code{\link{predict.glm}},\code{\link{PrecipitationOccurenceModel}},\code{\link{PrecipitationOccurenceMultiSiteModel}}
#' @export 
#' @method generate PrecipitationOccurenceModel
#' @S3method generate PrecipitationOccurenceModel
#' @aliases generate generate.PrecipitationOccurenceModel 
#' @rdname generate
#' @importFrom RGENERATE generate
#'
#' @references
#' D.S. Wilks (1998), Multisite Generalization of a Daily Stochastic Precipitation Generation Model, Journal of Hydrology, Volume 210, Issues 1-4, September 1998, Pages 178-191,
#' \url{http://www.sciencedirect.com/science/article/pii/S0022169498001863}
#' 
#' Muamaraldin Mhanna and Willy Bauwens (2011) A Stochastic Space-Time Model for the Generation of Daily Rainfall in the Gaza Strip, International Journal of Climatology, Volume 32, Issue 7, pages 1098-1112,
#' \url{http://dx.doi.org/10.1002/joc.2305}
#' 
#' 
#' 
#'  
#' @examples
#' 
#' library(RGENERATEPREC)
#' 
#' 
#' ## A function example can be found in the following script file: 
#' scriptfile <- system.file("example.generate.R",package="RGENERATEPREC")
#' ## The corrent file path is given by 'scriptfile' variable:
#' print(scriptfile)
#' ## To run the example file, launch the file with 'source' command (uncomment the following line)
#' #source(scriptfile)
#' 
#' ## ALTERNATIVELY you can run the following lines:
#' 
#' 
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

#' it <- station[2]
#' vect <- Tx_mes[,it]-Tn_mes[,it]
#' months <- factor(prec_mes$month)
#' 
#' #
#' ### Not Run!!! 
#' ###  Please uncomment the following lines to run them 
#'  
#' 
#' #model <- 
#' #PrecipitationOccurenceModel(x=prec_mes[,it],exogen=vect,
#' #monthly.factor=months,valmin=valmin)
#' #
#' #obs <- prec_mes[,it]>=valmin
#' #
#' #gen <- generate(model,exogen=vect,monthly.factor=months,n=length(months))
#' 
#' 
#' ### MultiSite Generation 
#' 
#' 
#' station <- station[1:2]
#' exogen <- Tx_mes[,station]-Tn_mes[,station]
#' 
#' months <- factor(prec_mes$month)
#' 
#' #
#' ### Not Run!!! 
#' ###  Please uncomment the following lines to run them 
#'  
#' #model_multisite <- 
#' #PrecipitationOccurenceMultiSiteModel(x=prec_mes[,station],
#' #exogen=exogen,origin=origin,multisite_type="wilks")
#' #
#' #
#' ## LOGIT-type Model 
#' #model_multisite_logit <- 
#' #PrecipitationOccurenceMultiSiteModel(x=prec_mes,exogen=exogen,
#' #origin=origin,multisite_type="logit",station=station)
#' #
#' #
#' #obs_multisite <- prec_mes[,station]>=valmin
#' #
#' #gen_multisite <- generate(model_multisite,exogen=exogen,origin=origin,end=end)
#' #
#' #gen_multisite_logit <- generate(model_multisite_logit,exogen=exogen,origin=origin,end=end)


generate.PrecipitationOccurenceModel <- function(x,newdata=NULL,previous=NULL,n=30,random=runif(n,min=0,max=1),exogen=NULL,monthly.factor=NULL,...) {
	
	p <- x$p
	
	if (p<1) previous <- NULL
	
	if (!is.null(exogen)) newdata <- as.data.frame(exogen)
	
	if (is.null(newdata) & is.null(monthly.factor)) {
		
		newdata <- x$predictor
		
	} else if (is.null(newdata)) {
		
		newdata <- as.data.frame(array(NA,c(length(monthly.factor),0)))
	}
	
	
	if (!is.null(monthly.factor)) newdata$month <- factor(monthly.factor)
	
	
	if (nrow(newdata)<n) {
		
		warning("Warning: n is reduced, insufficient numbers of predictors!")
		n <- nrow(newdata)
		
	}
	
	names_n <- names(newdata)
	newdata <- as.data.frame(newdata[1:n,])
	names(newdata) <- names_n
	
	if (is.null(previous)) {
		
		previous <- rnorm(x$p)>=0
		
	}
	
	
	out <- array(NA,n)
	
	
	for (i in 1:n) {
		

###		prob <- 1-predict(x,newdata=newdata[i,],previous=previous,type="response",endogenous=endogenous,...)
		prob <- 1-predict(x,newdata=newdata[i,],previous=previous,type="response",...) ## ec 20141208
		out[i] <- random[i]>=prob
		

		previous  <- c(out[i],previous[-p])

	}
		
	
	
	return(out)
} 


NULL
#'
#' 
#' 
#' @export 
#' @method generate CCGammaObjectListPerEachMonth
#' @S3method generate CCGammaObjectListPerEachMonth
#' @aliases generate generate.CCGammaObjectListPerEachMonth
#' @rdname generate
#' 

generate.CCGammaObjectListPerEachMonth <- function(x,...) {
	
	
	class(x) <- "list"
	out <- generate(x,...)
	
	return(out)
	
	
}





NULL
#'
#' 
#' 
#' @export 
#' @method generate PrecipitationOccurenceMultiSiteModel
#' @S3method generate PrecipitationOccurenceMultiSiteModel
#' @aliases generate generate.PrecipitationOccurenceMultiSiteModel
#' @rdname generate
#' 

generate.PrecipitationOccurenceMultiSiteModel <- function(x,exogen,n=10,origin="1961-1-1",end="1990-1-1",previous=NULL,monthly.factor=NULL,...) {
	
	
	out <- NULL 
	if (is.null(monthly.factor)) {
		
		dates <- as.Date(origin):as.Date(end)
		months <- adddate(as.data.frame(dates),origin=origin)$month
		
		n <- length(months)
		
		
	} else {
		months <- monthly.factor
		n <- length(months)
	}
	if (x$type=="wilks") {
		
		monthsf <- sprintf("month%02d",months)
		
		
		gen_wilks <- generate(x$ccgamma,FUN=rnorm,type="covariance",names=x$station,factor.series=monthsf)
	
		for (c in 1:ncol(gen_wilks)) {
			
			gen_wilks[,c] <- pnorm(gen_wilks[,c])
			
		} 
		
		
		
		
	 

	
		if (is.null(exogen)) {
			
			exogen <- lapply(X=x$station,FUN=function(x){ NULL })
			names(exogen) <- x$station
			
			
		}
		
		if (is.null(previous)) {
			
			previous <- lapply(X=x$station,FUN=function(x){ NULL })
			names(previous) <- x$station
			
		}
		out <- as.data.frame(array(NA,dim(gen_wilks)))
		names(out) <- names(gen_wilks)
		
		for (it in x$station) {
			
#			if (is.data.frame(exogen)) {
#				
#				exogen_loc <- exogen[,it]
#				
#			} else {
#				
#				exogen_loc <- exogen[[it]]
#			}
			
			if (is.data.frame(exogen)) {
				
				cols <- str_detect(names(exogen),it)
				exogen_loc <- exogen[,cols]
				
				
			} else if (is.list(exogen)) {
				
				exogen_loc <- exogen[[it]]
			} else {
				
				exogen_loc <- exogen
			}
			
			if (is.data.frame(previous)) {
				
				previous_loc <- previous[,it]
				
			} else {
				
				previous_loc <- previous[[it]]
			}
			
			### 
			###function(x,newdata=NULL,previous=NULL,n=30,random=runif(n,min=0,max=1),exogen=NULL,monthly.factor=NULL,...) {
			message(paste("Processing",it))
			out[,it] <- generate(x[[it]],previous=previous_loc,exogen=exogen_loc,monthly.factor=factor(months),random=gen_wilks[,it],n=n)
		
			 
			
			###
		} 
		
		
		
		
		
	}	else if (x$type=="logit") {
	
		if (is.null(exogen)) { 
		
			exogen <- as.data.frame(array(NA,c(n,0)))
			
		}	
		
		if (is.null(previous)) {
			
			
				
				previous <- as.data.frame(array(rnorm(x$p*x$K)>=0,c(x$p,x$K)))
				names(previous) <- x$station
				
				
			
		}	else {
			
			
			previous <- previous[,x$station] ## ec 20141204
		}
	
		out <- as.data.frame(array(NA,c(n,length(x$station))))
		names(out) <- x$station
		

		percs <- seq(from=0,to=100,by=5)
		npercs <- trunc(percs/100*n)
		
		
		

		for (ncnt in 1:n) {
			
			if (ncnt %in% npercs) {
				
				valprec <- percs[npercs==ncnt]
				message <- paste(sprintf("Processing: %0.2f",valprec),"%",sep="")
				message(message)
				
			}

		    out[ncnt,] <- unlist(lapply(X=x[x$station],FUN=generate,previous=previous,endogenous=x$station,exogen=exogen[ncnt,],monthly.factor=factor(months)[ncnt],n=1,...))

			previous[-1,] <- previous[-x$p,]
			previous[1,] <- out[ncnt,]
			
		}
		
		
		
		
####		out[,it] <- generate(x[[it]],previous=previous_loc,exogen=exogen,monthly.factor=factor(months),random=gen_wilks[,it],n=n)
		
		
		
		
	}
	###		out <- NULL 
	
	## TO DO 
	
	
	## TO GO ON ....
	
	return(out)
}


