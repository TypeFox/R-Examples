#' Iterator that returns the Cartesian product of the arguments.
#'
#' Constructs an iterator that is the Cartesian product of each of the arguments.
#'
#' Although they share the same end goal, \code{iproduct} can yield drastic
#' memory savings compared to \code{\link[base]{expand.grid}}.
#'
#' @importFrom iterators iter nextElem
#' @export
#' @param ... multiple arguments
#' @return iterator that iterates through each element from the Cartesian
#' product
#' 
#' @examples
#' it <- iproduct(x=1:3, y=4:5)
#' iterators::nextElem(it) # list(x=1, y=4)
#' iterators::nextElem(it) # list(x=1, y=5)
#' iterators::nextElem(it) # list(x=2, y=4)
#' iterators::nextElem(it) # list(x=2, y=5)
#' iterators::nextElem(it) # list(x=3, y=4)
#' iterators::nextElem(it) # list(x=3, y=5)
#'
#' # iproduct is a replacement for base::expand.grid()
#' # Large data.frames are not created unless the iterator is manually consumed
#' a <- 1:2
#' b <- 3:4
#' c <- 5:6
#' it2 <- iproduct(a=a, b=b, c=c)
#' df_iproduct <- do.call(rbind, as.list(it2))
#' df_iproduct <- data.frame(df_iproduct)
#'
#' # Compare df_iproduct with the results from base::expand.grid()
#' base::expand.grid(a=a, b=b, c=c)
#' 
iproduct <- function(...) {
  args_list <- list(...)
  if (length(args_list) == 0) {
    stop("At least one argument must be supplied.")
  }

  # Determines the frequency and number of replicates each element in args_list
  # must be repeated
  args_lengths <- sapply(args_list, length)
  freq <- cumprod(c(1, rev(args_lengths[-1])))
  freq <- rev(unname(freq))
  rep_times <- unname(prod(args_lengths) / freq / args_lengths)

  # For each element in args_list, an iterator is constructed that repeats each
  # of its values with the proper frequency and cadence.
  args_iters <- mapply(irep,
                       args_list,
                       times=rep_times,
                       each=freq,
                       SIMPLIFY=FALSE)

  # Finally, izip's the list of iterators
  it <- do.call(izip, args=args_iters)
}
