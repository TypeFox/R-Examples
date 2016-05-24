
STATUSLEVELS <- data.frame(
	short = c("success","no memory","diverges","max iterations"),
	long =  c(
      "success",
      "could not allocate enough memory",
      "divergence detected",
      "maximum number of iterations reached"
   ),
	stringsAsFactors=FALSE
)

# make ordered status vector
new_status <- function(n){
   statuslabels=STATUSLEVELS[,'long']
   ordered(statuslabels[n+1], levels=statuslabels )
  
}

#' Adjusted object
#' @name adjusted
#'
#' @seealso \code{\link{adjust}}
#'
#' @section Details: 
#' A \code{adjusted} object contains the adjusted vector as well as some information on how
#' the adjustment was achieved. In particular, it contains the following slots (to be accessed with
#' the dollar operator):
#'    \itemize{
#'       \item{ \code{\$x}: the adjusted vector.}
#'       \item{ \code{\$accuracy}: Maximum deviance of \code{\$x} from the constraints (see \code{\link{adjust}} for details).}
#'       \item{ \code{\$objective} : Square root of objective function \eqn{\sum_i(x_i-x^0_i)^2w_i}.}
#'       \item{ \code{\$duration}: \code{proc_time} object showing time it took to run the adjustment. (See \code{proc.time}).}
#'       \item{ \code{\$niter}: Number of iterations.}
#'       \item{ \code{\$status}: A \code{character} string stating whether the adjustment was successful, 
#'                aborted, or if the maximum number of iterations was reached before convergence. }
#'       \item{ \code{\$method}: \code{'sparse'} or \code{'dense'}.}
#'    }
#'
#'
#'
{}


#' @method print adjusted
#' @param x an object of class \code{adjusted}
#' @param maxprint max number of output values to print
#' @param ... parameters to pass to other methds
#' 
#' @rdname adjusted
#' @export
print.adjusted <- function(x, maxprint = 10, ...){
    cat("Object of class 'adjusted'\n")
    cat(sprintf("  Status    : %s (using '%s' method)\n", x$status, x$method))
    cat(sprintf("  Accuracy  : %g\n", x$accuracy))
    cat(sprintf("  Objective : %g\n", x$objective))
    cat(sprintf("  Iterations: %d\n", x$niter))
    cat(sprintf("  Timing (s): %g\n", x$duration['elapsed']))
    tr = ":";
    if (length(x$x) > 10) tr = sprintf(" (truncated at %d):",maxprint)
    cat(paste("Solution",tr,"\n",sep=""))
    print(x$x[1:min(length(x$x),maxprint)])
}

# create 'adjusted' object. Input is a solution vector 
# with attributes, returned by "R_sc_solve_spa" or "R_dc_solve_spa"
new_adjusted <- function(x, duration, method, objective, varnames=NULL){
   acc = attr(x,"accuracy")
   nit = attr(x,"niter")
   status = new_status(attr(x,"status"))
   attr(x,"accuracy") <- NULL
   attr(x,"niter")    <- NULL
   attr(x,"status")   <- NULL
   names(x) <- varnames
   structure(
      list(x = x, accuracy = acc, objective=objective, method=method, niter = nit, duration=duration, status=status ),
      class = "adjusted"
   )   
}

