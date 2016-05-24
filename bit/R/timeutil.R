# timing utilities for ff and bit
# (c) 2012 Jens Oehlschlägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2012-05-28

#! \name{repeat.time}
#! \alias{repeat.time}
#! \title{
#! Adaptive timer
#! }
#! \description{
#! Repeats timing expr until minSec is reached
#! }
#! \usage{
#! repeat.time(expr, gcFirst = TRUE, minSec = 0.5, envir=parent.frame())
#! }
#! \arguments{
#!   \item{expr}{Valid \R expression to be timed.}
#!   \item{gcFirst}{Logical - should a garbage collection be performed
#!     immediately before the timing?  Default is \code{TRUE}.}
#!   \item{minSec}{number of seconds to repeat at least}
#!   \item{envir}{the environment in which to evaluate \code{expr} (by default the calling frame)}
#! }
#! \value{
#!   A object of class \code{"proc_time"}: see
#!   \code{\link{proc.time}} for details.
#! }
#! \seealso{
#!   \code{\link{system.time}}
#! }
#! \author{
#! Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \examples{
#!   system.time(1+1)
#!   repeat.time(1+1)
#!   system.time(sort(runif(1e6)))
#!   repeat.time(sort(runif(1e6)))
#! }
#! \keyword{utilities}

repeat.time <- 
function (expr, gcFirst = TRUE, minSec=0.5, envir=parent.frame()) 
{
	ppt <- function(y) {
		if (!is.na(y[4L])) 
			y[1L] <- y[1L] + y[4L]
		if (!is.na(y[5L])) 
			y[2L] <- y[2L] + y[5L]
		y[1L:3L]
	}
	if (!exists("proc.time")) 
		return(rep(NA_real_, 5L))
	if (gcFirst) 
		gc(FALSE)
	time <- proc.time()
	on.exit(cat("Timing stopped at:", ppt(proc.time() - time), 
		"\n"))
	r <- 0L
	while((proc.time()[3]-time[3]) < minSec){
		r <- r + 1L
		eval(substitute(expr), envir=envir)
	}	
	new.time <- proc.time()
	on.exit()
	structure((new.time - time)/r, class = "proc_time")
}	
