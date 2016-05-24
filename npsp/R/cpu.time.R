#--------------------------------------------------------------------
#   cpu.time.R (npsp package)
#--------------------------------------------------------------------
#   .cpu.time.ini()
#   cpu.time(..., reset, total, last, flush)
#--------------------------------------------------------------------
# CPU time utilities
#   Call cpu.time(restart = TRUE) where you want to start counting.
#   Call cpu.time() to print/get total and/or partial (since the last call 
#   to this function) real and CPU times .
#--------------------------------------------------------------------

#--------------------------------------------------------------------
#' @rdname npsp-internals
#' @keywords internal
#' @export
.cpu.time.ini <- function() {
    time.ini <- structure(rep(0, 5), .Names = c("user.self", "sys.self", "elapsed",  
        "user.child", "sys.child"), class = "proc_time")# proc.time()
    time.last <- time.ini
    function(..., reset = FALSE, total = TRUE, last = TRUE, flush = FALSE) {
        res <- list(time = proc.time())
        if (reset) {
            time.ini <<- res$time
            time.last <<- time.ini        
            res$last <- res$total <- 0
            if (total | last) cat("CPU time has been initialized.\n")
        } else {        
            res$last <- res$time - time.last
            res$total <- res$time - time.ini
            if (last) {
                cat("Time of last operation:", ..., "\n")
                print(res$last)
            }    
            if (total) {
                cat("Total time:\n")
                print(res$total)
            }
            if (flush) flush.console()
            time.last <<- res$time
        }    
        return(invisible(res))
    }
#--------------------------------------------------------------------
}

#--------------------------------------------------------------------
#' Total and partial CPU time used
#' 
#' Returns and (optionally) prints the total and/or partial (since the last call to this function) 
#' real and CPU times.
#' @param ... objects (describing the last operation) to be printed (using \code{\link{cat}}), 
#' if \code{last == TRUE}.
#' @param reset logical; if \code{TRUE}, time counters are initialized. 
#' @param total logical; if \code{TRUE}, the total time used is printed. 
#' @param last logical; if \code{TRUE}, the partial time used is printed. 
#' @param flush logical; if \code{TRUE}, \code{\link{flush.console}} is called.
#' @return Invisibly returns a list with  the following 3 components 
#' (objects of class \code{"proc_time"}):
#' \item{time}{user, system, and total elapsed times for the currently running R process 
#' (result of a call to \code{\link{proc.time}}). }
#' \item{last, total}{differences between the corresponding \code{\link{proc.time}} calls.}
#' @seealso
#' \code{\link{proc.time}}, \code{\link{system.time}}, \code{\link{flush.console}}.
#' @examples
#' ## Not run:
#'  
#'  cpu.time(reset=TRUE)
#'  res <- median(runif(100000))
#'  cpu.time('\nSample median of', 100000, 'values =', res)
#'  res <- median(runif(1000))
#'  cpu.time('\nSample median of', 1000, 'values =', res)
#' 
#' ## End(Not run) 
#' @export
cpu.time <- .cpu.time.ini()


