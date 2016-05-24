#' Add Memoization to a Function
#'
#' @param fcn the function to be memoized
#' @param verbose show debugging information
#' 
#' @examples 
#' runifCached <- addMemoization(runif)
#'
#' @return memoized function
#' @export
addMemoization <- function(fcn, verbose=FALSE) {
  if(verbose) {
    cat("In modified addMemoization()", sep = "\n")
  }
  
  if (!is.function(fcn)) {
    stop("Argument 'fcn' is not a function");
  }

  # Already memoized?
  if (inherits(fcn, "MemoizedFunction")) {
    return(fcn)
  }

  res <- function(...) {
    args <- list(fcn, ...)
    do.call("memoizedCall", args=args)
  }
  class(res) <- c("MemoizedFunction", class(res))

  res
} 
