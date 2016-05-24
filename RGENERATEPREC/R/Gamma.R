#
#
NULL
#'
#' This function extends \code{\link[RMAWGEN]{continuity_ratio}} and adds the corresponding gaussian correlation matrix for no-precipitation occurence. 
#'
#' @param data data frame or 'zoo'  R object containing daily precipitation time series for several gauges (one gauge time series per column). See \code{\link[RMAWGEN]{continuity_ratio}}.
#' @param lag numeric lag (expressed as number of days) used for computation for "cross" continuity ratio and joint probability of prercipitation (no)occurence. See \code{\link[RMAWGEN]{continuity_ratio}}.
#' @param p positive integer parameter. Default is \code{NA}, otherwise, \code{lag} is calculated as the vector \code{0:p}. 
#' @param valmin threshold precipitation value [mm] for wet/dry day indicator. 
#' If precipitation is lower than \code{valmin}, day is considered dry. Default is 0.5 mm. See \code{\link[RMAWGEN]{continuity_ratio}}.
#' @param p0_v1 vector for marginal probablities, see \code{\link{omega}} and \code{\link{omega_inv}}. 
#' @param nearPD see \code{\link{omega_inv}}. Default is \code{(lag==0)}.
#' @param interval,tolerance see \code{\link{omega_inv}}
#' @param only.matrix logical value. If \code{TRUE} the function returns only the gaussian correlaton matrix. Deafaul is \code{FALSE}.
#' @param return.value string. If it is not either \code{NULL} (Default) and \code{NA}, function returns only the argument indicated by this argument.
#' @param null.gcorrelation numerical value \code{nooccurence_gcorrelation}  under which is considered to be 0. 
#' @param sample character string indicated if function must be calculated differently for subset of the year, e.g. monthly. Admitted values are \code{NULL} (Default), \code{"all"} or \code{"monthly"}. 
#' @param origin character string (yyyy-dd-mm) indicated the date of the first row of \code{"data"}. It is used if \code{data} and \code{sample} are not \code{NULL}. 
#' @param ... additional agruments of \code{\link{omega_inv}} or \code{\link{CCGamma}}

#' @author Emanuele Cordano
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
#' 
#' @return  An object which is a list containing the following fields: 
#'
#' \code{continuity_ratio} : \code{lag}-day lagged  continuity ratio, as returned by \code{\link[RMAWGEN]{continuity_ratio}}; 
#'
#' \code{occurence}  : joint probability of \code{lag}-day lagged precipitation occurence, as returned by \code{\link[RMAWGEN]{continuity_ratio}}; 
#' 
#' \code{nooccurence} : joint probability of \code{lag}-day lagged no precipitation occurence, as returned by \code{\link[RMAWGEN]{continuity_ratio}};
#' 
#' \code{lag} : number of days lagged between the two compared events (see argument \code{lag});
#' 
#' \code{p0_v1} : vector of marginal probability of no precipitation occurence. If \code{lag}
#' is 0, it corresponds to the diagonal of  \code{nooccurence} matrix (see argument \code{p0_v1});
#' 
#' \code{nooccurence_gcorrelation} corresponding gaussian correlation for no precipitation occurence obtained by applying \code{\link{omega_inv}} to \code{nooccurence},
#' 
# ##### \code{TransintionMatrixMCFirstOrder} List of transition matrices (with condinitioned probability values) of the first-order Markov Chain.  See \code{\link{TransitionMatrixMCFirstOrder}}.
#' 
#' If the argument \code{only.matrix} is \code{TRUE}, only \code{nooccurence_gcorrelation} is returned as a matrix. 
#' In case the argument \code{lag} is a vector wirh length more than one, the function returns a list of the above-cited return object for each value of the vector \code{lag}.
#' 
#' 
#' @note This functon is useful to generate the serial cross-correlation matrices for no precipitation occurrence for Yule-Walker Equations. In case \code{lag} is a vactor, \code{nearPD} must be a vector of the same size, 
#' default is \code{(lag==0)}.
#'  
#' See the R code for major details
#' 
#' @seealso \code{\link[RMAWGEN]{continuity_ratio}},\code{\link{omega_inv}},\code{\link{omega}},\code{\link{CCGammaToBlockmatrix}}
#' @importFrom RMAWGEN adddate
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
#' CCGamma <- CCGamma(data=prec_mes,lag=0,tolerance=0.001,only.matrix=FALSE)
#' ## Not Run in the examples, uncomment to run the following line
#' ## CCGamma <- CCGamma(data=prec_mes,lag=0:2,tolerance=0.001,only.matrix=FALSE)
#' 
#' ## Not Run in the examples, uncomment to run the following line
#' ## CCGamma_monthly <- CCGamma(data=prec_mes,lag=0,tolerance=0.001,only.matrix=FALSE,
#' #                     sample="monthly",origin=origin)
#' 
#' 
#' @export


# INSERIRE EXAMPLES!!!  nearPD=(lag==0)

CCGamma <- function(data,lag=0,p0_v1=NULL,p=NA,valmin=0.5,nearPD=(lag>=0),interval=c(-1,1),tolerance=.Machine$double.eps,only.matrix=FALSE,return.value=NULL,null.gcorrelation=1e-5,sample=NULL,origin="1961-1-1",...)
	
			{
				
				out <- NULL
				if (is.null(sample)) sample <- "all"

				if (sample=="monthly") {
						
					out <- list()
					data <- adddate(data,origin=origin)
					month <- data$month
					months <- unique(data$month)
					
					data <- data[,!(names(data) %in% c("year","month","day"))]
					for (m in months) {
						
				#		mprev <- m-1
				#		mnext <- m+1
				#		mprev[mprev==0] <- 12
				#		mnext[mnext==13] <- 1
						
						cond <- month==m
					#	cond <- ((data$month==mprev) & (data$day>=25)) | cond
					#	cond <- ((data$month==mnext) & (data$day<=8)) | cond	   
					#	data <- data[cond,]
					#	data <- data[,names]
			
						out[[m]] <- CCGamma(data=data[cond,],
								lag=lag,p0_v1=p0_v1,p=p,valmin=valmin,nearPD=nearPD,interval=interval,tolerance=tolerance,only.matrix=only.matrix,return.value=return.value,null.gcorrelation=null.gcorrelation,sample="all",origin=NULL,...)
						
						
						
						
						
					}
					
					names(out) <- months 
					
					class(out) <- "CCGammaObjectListPerEachMonth"
					
					return(out)
				}
				
				
				
				
				out <- NULL 
				text <- "nooccurence"
				
				if ((length(p0_v1)==0)) {
					
					temp <- continuity_ratio(data,lag=0,valmin=valmin)
					p0_v1 <- diag(temp[[text]])
					
				} 
				
				if (!is.na(p)) {
					
					lag <- 0:p
					nearPD <- (lag==0)
				}
				if (length(lag)==1) {
					out <- continuity_ratio(data,lag=lag,valmin=valmin)
					
	
					### NOOCCURENCE_CORRELATION 
					out$lag <- lag
					out$p0_v1 <- p0_v1
					message("lag")
					message(lag)
					out$nooccurence_correlation <- out[[text]]
					
					#### CORRECTION HERE FOR p0v
					
					
					
					#####
					for (r in 1:nrow(out$nooccurence_correlation)) {
						for (c in 1:ncol(out$nooccurence_correlation)) {
							
							p0v <- diag(out[[text]])
							out$nooccurence_correlation[r,c] <- (out[[text]][r,c]-p0_v1[r]*p0_v1[c])/(p0_v1[r]*p0_v1[c]*(1-p0_v1[r])*(1-p0_v1[c]))^0.5
							
						}
						
					}
					
					
					### NOOCCURENCE_CORRELATION
				    
					
					out$nooccurence_gcorrelation <- omega_inv(p0=out[[text]],p0_v1=p0_v1,correlation=NA,only.value=TRUE,interval=interval,tolerance=tolerance,nearPD=nearPD,...)
					out$nooccurence_gcorrelation[abs(out$nooccurence_gcorrelation)<null.gcorrelation] <- 0 ## added on 2014 03 18
			
					### out$TransintionMatrixMCFirstOrder
					
					####out$TransintionMatrixMCFirstOrder <- TransitionMatrixMCFirstOrder(as.data.frame(data>valmin),rc.names=c("dry","wet"))
		           
		
					if (only.matrix) {
						out <- out$nooccurence_gcorrelation
					} else {
						
						class(out) <- "CCGammaObject"
					}
				###	if (!is.null(return.value)) return.value <- NA
				
				
					if (!is.null(return.value)) { out <- out[[return.value]] }
				} else if (length(lag)>1) {
				
					out <- as.list(array(NA,length(lag)))
					names(out) <- lag
				
					for (i in 1:length(lag)) {
						
				
						out[[i]] <- CCGamma(data=data,lag=lag[i],nearPD=nearPD[i],p0_v1=p0_v1,valmin=valmin,interval=interval,tolerance=tolerance,only.matrix=only.matrix,return.value=return.value,...)
					} 
					
					
					class(out) <- "CCGammaObjectList"
				}
				
				
				return(out)
			
					
			}