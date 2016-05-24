#' All information used to calculate Kendall's distance.
#' 
#' Performs each column-wise comparison on a matrix of sequences.  A 0 value
#' denotes that there is an increase between the two columns, 1 a decrease, and
#' NA indicates that the column values are identical in the row.
#' 
#' 
#' @param r Matrix of sequences.
#' @param inds Possibly efficiency increase when doing repeated calculations,
#' currently not used.
#' @return Matrix of 0s, 1s, and NAs representing pairwise comparisons of
#' vector values.
#' @author Erik Gregory
#' @references http://en.wikipedia.org/wiki/Kendall_tau_distance
#' @keywords Kendall Distance

KendallInfo <-
function(r, inds = NULL) {
  if (is.null(inds)) {
    inds <- combn(ncol(r), 2)
  }
  if (!is.matrix(r)) {
    r <- as.matrix(r)
    attr(r, "dimnames") <- NULL
  }
  infos <- r[, inds[1, ]] - r[, inds[2, ]]
  decr <- which(infos > 0)
  incr <- which(infos < 0)
  indet <- which(infos == 0)
  infos[decr] <- TRUE
  infos[incr] <- FALSE
  infos[indet] <- NA
  return(infos)
}
