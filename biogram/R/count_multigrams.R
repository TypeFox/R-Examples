#' Detect and count multiple n-grams in sequences
#'
#' A convinient wrapper around \code{\link{count_ngrams}} for counting multiple
#' values of \code{n} and \code{d}.
#'
#' @inheritParams count_ngrams
#' @param ns \code{numeric} vector of n-grams' sizes. See Details.
#' @param ds \code{list} of distances between elements of n-grams. Each element of the list
#' is a vector used as distance for the respective n-gram size given by the \code{ns}
#' parameter.
#' @return a \code{integer} matrix with named columns. The naming conventions are the same
#' as in \code{\link{count_ngrams}}.
#' @details \code{ns} vector and \code{ds} vector must have equal length. Elements of 
#' \code{ds} vector are used as equivalents of \code{d} parameter for respective values 
#' of \code{ns}. For example, if \code{ns} is \code{c(4, 4, 4)}, the \code{ds} must be a list of 
#' length 3. Each element of the \code{ds} list must have length 3 or 1, as appropriate
#' for a \code{d} parameter in \code{count_ngrams} function.
#' @export
#' @examples 
#' seqs <- matrix(sample(1L:4, 600, replace = TRUE), ncol = 50)
#' count_multigrams(c(3, 1), list(c(1, 0), 0), seqs, 1L:4, pos = TRUE)
#' #if ds parameter is not present, n-grams are calculated for distance 0
#' count_multigrams(c(3, 1), seq = seqs, u = 1L:4)
#' 
#' #calculate three times n-gram with the same length, but different distances between
#' #elements
#' count_multigrams(c(4, 4, 4), list(c(2, 0, 1), c(2, 1, 0), c(0, 1, 2)), 
#'                  seqs, 1L:4, pos = TRUE)

count_multigrams <- function(ns, ds = rep(0, length(ns)), seq, u, pos = FALSE, scale = FALSE, threshold = 0) {
  if(length(ns) != length(ds))
    stop("'ns' vector must have equal length to 'ds' list.")
  
  n_loops <- length(ns)
  do.call(cbind, lapply(1L:n_loops, function(current_loop) {
      count_ngrams(seq, ns[current_loop], u, ds[[current_loop]], 
                   pos, scale, threshold)
  }))
}


