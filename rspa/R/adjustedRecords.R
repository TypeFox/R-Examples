#' Adjusted records
#' @name adjustedRecords
#' @seealso \code{\link{adjustRecords}}
#' @section Details:
#' The \code{adjustedRecords} object contains adjusted data as well as some information 
#' on the adjusting process. In particular:
#' \itemize{
#'		\item{\code{\$adjusted}: \code{data.frame} similar to \code{dat} with adjusted values.}
#'	
#'    \item{\code{$status}: \code{data.frame} with the same number of rows as \code{dat}. Each row
#'			stores the information of one \code{\link{adjusted}} object.}
#' }
#'
#' When printed, only the first 10 rows of \code{cbind(adjusted, status)} are shown. Use \code{summary}
#' for a quick overview. The \code{plot} method shows kernel density estimates of the accuracy and
#' objective functions. To avoid densities at values below 0, the accuracy densities are evaluated
#' under a sqrt-transform and transformed back before plotting. For the objective function values
#' a log-transform is used.
#' 
{}


#' @method print adjustedRecords
#' @param x object of class \code{adjustedRecords}
#' @param ... additional parameters to pass to other methods
#' @rdname adjustedRecords
#' @export
print.adjustedRecords <- function(x,...){
   I <- 1:min(10,nrow(x$adjusted))
   cat("Object of class 'adjustedRecords'\n")
   print(
      cbind(x$adjusted[I,],x$status[I,])
   )
   if ( nrow(x$adjusted)>10 ) cat("print truncated...\n")
}


#' @method summary adjustedRecords
#' @param object object of class \code{adjustedRecords}
#' @rdname adjustedRecords
#' @export
summary.adjustedRecords <- function(object,...){
	cat("Object of class 'adjustedRecords'\n")

   lv <- levels(object$status$status)
	nsuccess <- sum(object$status$status == lv[1],na.rm=TRUE)
	p <- list()
	p[[1]] <- sprintf(" Records : %d\n",nrow(object$adjusted))
	p[[2]] <- sprintf(" Adjusted: %d (%d converged) \n",sum(!is.na(object$status$status)),nsuccess)
	p[[3]] <- sprintf(" duration: %gs (total)\n",sum(object$status$elapsed)) 
	p[[4]] <- ""
	p[[5]] <- "Summary of adjusted records:\n"
	d <- lapply(p,cat)
	iadj <- !is.na(object$status$status)
   print(summary(object$status[iadj,c('objective','accuracy')]))
}


#'
#' @method plot adjustedRecords
#' @rdname adjustedRecords
#' @export
plot.adjustedRecords <- function(x,...){

   if (nrow(x$adjusted) <= 1 ){
		cat("Nothing to plot...\n")
		return();
	}

	par(mfrow=c(2,1),mar=c(2,4,4,1))
	lwd = '2'


   
   ia <- !is.na(x$status$status) & x$status$accuracy > 0
	a <- x$status$accuracy[ia]
	if ( length(a) >=2 ){
		d <- density(sqrt(a))
		plot(d$x*d$x, d$y + .Machine$double.eps,
			main= sprintf("Accuracy (%d of %d adjusted records positive)",length(a),sum(!is.na(x$status$status))),
			ylab='density',
			xlab='',
			type='l',
			lwd=lwd
		)
		rug(a,col="blue")
	} else {
		plot.new()
		text(0.5,0.5,"too few points to plot density")
	}

	a <- x$status$objective[x$status$objective > 0]
	if ( length(a) >= 2 ){
		d <- density(log(a))	
		plot(exp(d$x), d$y + .Machine$double.eps,
			main= sprintf("Objective function (%d of %d records adjusted)",length(a),sum(!is.na(x$status$status))),
			ylab='density',
			xlab='',
			type='l',
			lwd=lwd,
			log='x'
		)
		rug(a,col="blue")
	} else {
		
		plot.new()
		text(0.5,0.5,"too few points to plot density")
	
	}

}



# internal method for statusblock of adjustedRecords
`%++%` <- function(x=NULL,y){
   if ( is.null(x) ) return(y)

   x$accuracy <- pmax(x$accuracy, y$accuracy)
   x$objective <- sqrt(x$objective*x$objective + y$objective*y$objective)
   x$status <- pmax(x$status,y$status,na.rm=TRUE)
   x$niter    <- x$niter + y$niter
   x[5:7]     <- x[5:7] + y[5:7] # durations
   x

}




