NULL 
#' 
#'\code{arch.test} function for \code{varest2} object
#'
#' @param object a \code{varest2} object
#' @param interval string or subset interval of time (e.g. days) or length of this subset interval to which the ARCH test is applied (see Note). Default is \code{NULL}.  
#' @param overlap number of time instants (e.g. days) which are overlapped on two different subsequent intervals. Default is 20. It is used only if \code{interval} has length 1. 
#' @param list.output logical value. If \code{TRUE} the function returns a list of the test results of each interval. It is used if \code{interval} is not \code{NULL}. Default is \code{FALSE}. 
#' @param ...   further arguments  for \code{\link{arch.test}}
#' 
#' @export 
#'
#' @details This function is a wrapper of \code{\link{arch.test}}. It can compute the test also for some subsets (intervals) of the time-series or for all the time-series divided in overlapping intervals. The intervals considered for the ARCH test are defined with the argument \code{interval}. If \code{interval} is an integer number instead of a vector, it indicates the length of the intervals in which the time-series is split. If \code{interval} is set to \code{NULL}, the test is done on the comprehensive residual time-series without splitting.   
#' 
#' @seealso \code{\link{arch.test}}
#' @return One object or a list of objects with class attribute \code{varcheck} as reported in \code{\link{arch.test}}




arch_test <- function(object,interval=NULL,overlap=20,list.output=FALSE,...) {
	
	temp <- object@VAR
	
	interval_min <- 100
	interval_sip <- overlap
	res <- residuals(object)
	nrow <- nrow(res)
	out <- NULL
	
	if (is.null(interval)) interval <- 1:nrow
	
	if (length(interval)==1) {
		
		if (is.na(interval)) interval <- interval_min
		interval[interval>nrow] <- nrow
		interval[interval<(overlap+1)] <- overlap+1
		interval[interval<interval_min] <- interval_min


		start <- 1
		end <- interval
		cnt <- 1
		while (end[cnt] < nrow) {
			cnt <- cnt+1
			start[cnt] <- end[cnt-1]-interval_sip
			end[cnt] <- start[cnt]+interval-1
		}
		
		if (list.output) out <- list()
		for (i in 1:length(start)) {
	
			
			at <- arch_test(object,interval=start[i]:end[i],...)
			
			if (is.null(out)) out <- at
			if (list.output) { 
				name <- paste(start[i],end[i],sep=":")
				out[[name]] <- at 
			} else if (at$arch.mul$p.value<out$arch.mul$p.value) {
				out <- at
			}
			
		}	
		
	
		
		
		
	} else {
		
		
		interval <- interval[interval %in% 1:nrow(res)]
		var_in_arch_test <- VAR_mod(res[interval,],p=0)
		out <- arch.test(var_in_arch_test,...)
		
	}	
#	if (class(object)=="GPCAvarest2") {
#		if (length(object@GPCA_residuals)>0) {
#			temp <- VAR_mod(object@GPCA_residuals$final_results,p=0)
#		}
#	}
	
	
#	### class(temp) <- "varest"
	return(out)
	
}


