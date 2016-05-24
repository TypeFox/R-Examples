NULL
#'
#' Transition Matrix for a 1st order Markov Chain 
#'  
#' @param data precipitation occurrence data 
#' @param rc.names names of the states. 
#' @export
#' 
#' @references \url{http://stats.stackexchange.com/questions/36099/estimating-markov-transition-probabilities-from-sequence-data?rq=1}
#' 
#' @note Function rturn \code{NULL} if \code{data} contains \code{NA}s values.
#' @seealso \code{\link{CCGamma}}
#' 
#'
#' @examples 
#' 
#'  ## # Not Run in the examples, uncomment to run the following lines
#' 
#' #library(RMAWGEN)
#' #data(trentino)
#' #
#' #year_min <- 1961
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
#' #valmin <- 0.5
#' #
#' #mt <- TransitionMatrixMCFirstOrder(data=prec_mes>valmin,rc.names=c("dry","wet"))
#' #
#' #CCGamma <- CCGamma(data=prec_mes,lag=0,tolerance=0.001,valmin=valmin,only.matrix=FALSE)
#' #
#' #i <- 1
#' #pd <- CCGamma$p0_v1[i]
#' #pdv <- mt[[i]]
#' #pdcalc <- pd*pdv["dry","dry"]+(1-pd)*pdv["wet","dry"]
#' #
#' #pw <- 1-pd
#' #pwcalc <- pd*pdv["dry","wet"]+(1-pd)*pdv["wet","wet"]
#' #
#' ## verify
#' #pd-pdcalc
#' #pw-pwcalc
#' #

#x <- c(1,2,1,1,3,4,4,1,2,4,1,4,3,4,4,4,3,1,3,2,3,3,3,4,2,2,3)
#p <- matrix(nrow = 4, ncol = 4, 0)
#for (t in 1:(length(x) - 1)) p[x[t], x[t + 1]] <- p[x[t], x[t + 1]] + 1
#for (i in 1:4) p[i, ] <- p[i, ] / sum(p[i, ])

TransitionMatrixMCFirstOrder <- function(data,rc.names=NULL) {
	
	
	if (length(which(is.na(as.matrix(data))))>0) { 
	
		warning("TransitionMatrixMCFirstOrder: argument data contains NAs, function is not calculated!! ")
		return(NULL)
		
	}
	out <- list()	
	
	
	for (c in 1:ncol(data)) {
		
		out[[c]] <- NULL
		
		x <- as.vector(data[,c])
		x[is.na(x)] <- "NA"
		x <- factor(x)
		l <- levels(x)
		
		
		p <- array(0,c(length(l),length(l)))
		
		if (is.null(rc.names)) rc.names <- l 
		
		for (t in 1:(length(x) - 1)) p[x[t], x[t + 1]] <- p[x[t], x[t + 1]] + 1
		for (i in 1:length(l)) p[i, ] <- p[i, ] / sum(p[i, ])
		
		rownames(p) <- rc.names
		colnames(p) <- rc.names
		
		out[[c]] <- p
		
		
		
		
		
	}
	
	names(out) <- names(data)
	return(out)
}