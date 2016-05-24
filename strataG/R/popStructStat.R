#' @name popStructStat
#' @title Population structure statistics
#' @description Population structure statistics
#'
#' @param g a \linkS4class{gtypes} object.
#' @param nrep number specifying number of permutation replicates to use for 
#'   permutation test.
#' @param strata.mat an optional matrix of permuted stratifications. See Notes 
#'  for more details. Ignored if \code{nrep} is not \code{NULL}. 
#' @param keep.null logical. Keep the null distribution from the 
#'   permutation test?
#' @param prime.type type of G'st to calculate. Can be "nei" or "hedrick".
#' @param model,gamma,pairwise.deletion parameters passed to 
#'   \code{\link[ape]{dist.dna}}. Note that defaults for these arguments 
#'   (in particular \code{model}) are the same as in \code{dist.dna}.
#' @param ... optional arguments passed to or from other functions.
#'
#' @return A list with three elements:
#'   \tabular{ll}{
#'     \code{stat.name} \tab the name of the statistic.\cr
#'     \code{result} \tab a vector of the statistic estimate and the p-value, 
#'       if replicates were conducted.\cr
#'     \code{null.dist} \tab a vector of the null distribution from the 
#'       permutations.\cr
#'    }
#' 
#' @note If \code{strata.mat} is provided, it must be a numeric matrix of 
#'   integers from \code{0} to \code{k - 1}, where \code{k} is the number of 
#'   strata. Each column is a separate permutation and the first column is assumed 
#'   to represent the original stratification. If not provided \code{= NULL}, 
#'   stratification is taken from \code{g}. This argument is primarily used 
#'   internally by \code{\link{popStructTest}}.
#'
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @useDynLib strataG
#' @importFrom Rcpp sourceCpp
#'
NULL

# create matrix of permuted stratifications
.permStrata <- function(g, nrep = NULL) {
  st <- strata(g)
  if(any(is.na(st))) stop("cannot permute with unstratified samples.")
  strata.num <- cbind(as.numeric(st) - 1)
  if(is.null(nrep)) nrep <- 0
  if(nrep < 1) return(cbind(strata.num))
  null.perm <- sapply(1:nrep, function(i) sample(strata.num))
  return(cbind(strata.num, null.perm))
}

# check that matrix of permuted stratifications is of right format
.checkStrataMat <- function(strata.mat, g, nrep) {
  if(is.null(strata.mat)) {
    strata.mat <- .permStrata(g, nrep)
  } 
  if(!(is.matrix(strata.mat) & is.numeric(strata.mat))) {
    stop("'strata.mat' must be a numeric matrix")
  }
  if(nrow(strata.mat) != nInd(g)) {
    stop("'nrow(strata.mat)' is not equal to 'nInd(g)'")
  }
  zero.ref <- apply(strata.mat, 2, function(x) any(x == 0))
  if(!all(zero.ref)) {
    stop("'all columns in 'strata.mat' must have strata designations starting with 0")
  }
  unstrat <- apply(strata.mat, 2, function(x) any(is.na(x)))
  if(any(unstrat)) stop("there are unstratified samples in 'strata.mat'")
  return(strata.mat)
}

# create returned list from result vector
.formatResult <- function(result, stat.name, keep.null) {
  if(length(result) == 1) keep.null <- FALSE
  p.val <- if(length(result) == 1) NA else mean(result >= result[1], na.rm = TRUE) 
  if(is.nan(p.val)) p.val <- NA
  return(list(
    stat.name = stat.name, 
    result = c(estimate = result[1], p.val = p.val),
    null.dist = if(keep.null) result[-1] else NULL
  ))
}