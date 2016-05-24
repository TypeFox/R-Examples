#' @include spikeSlabGAM.R
{}

#' Convert samples from a model fitted with \code{spikeSlabGAM} into a
#' \code{bugs}-object
#'
#' Use \code{plot} and \code{print} on the returned
#' \code{\link[R2WinBUGS]{bugs}} object for convergence diagnostics. Lack of
#' convergence for \eqn{\alpha} or \eqn{\xi} does not necessarily mean that the
#' sampler has not converged for \eqn{\beta}.
#'
#' @param m a model fitted with \code{spikeSlabGAM}
#' @param rm which parameter blocks to omit from the conversion. Defaults to
#'   omission of \eqn{\alpha, \xi}, and \eqn{\gamma}. Set to NULL to keep all
#'   MCMC samples.
#' @return a \code{\link[R2WinBUGS]{bugs}} object which has convenience
#'   functions for assessing convergence based on parallel MCMC runs
#' @importFrom R2WinBUGS as.bugs.array
#' @export
#' @author Fabian Scheipl
#' @examples
#' #see help for spikeSlabGAM
ssGAM2Bugs <- function(m, rm = c("alpha", "ksi","gamma")) {
  stopifnot(class(m)=="spikeSlabGAM")

  list2array <- function(x) {
    a <- array(0, dim = c(nrow(x[[1]]), length(x), ncol(x[[1]])),
      dimnames = list(n.keep = NULL, n.chains = 1:length(x), colnames(x[[1]])))
    for(i in 1:length(x)) a[, i,] <- x[[i]]
    return(a)
  }
  smpls <- if(is.null(rm)) {
    m$samples
  } else {
    m$samples[!(names(m$samples) %in% rm)]
  }

  arrs <- lapply(smpls, list2array)

  #replace .digit at the end of a parameter name by [digit] so bugs recognizes
  # the grouping of the parameters:
  bugsnms <- gsub("\\[\\.", "\\[", gsub("(\\.[0-9]+$)","\\[\\1\\]",
    unlist(sapply(arrs, function(x) dimnames(x)[[3]]))))

  arr <-  array(0, dim = c(dim(arrs[[1]])[1:2],
    sum(sapply(arrs, function(x) dim(x)[3]))),
    dimnames = list(n.keep = NULL, n.chains = 1:m$mcmc$nChains, bugsnms))
  start <- 1
  for(i in 1:length(arrs)) {
    arr[,, start:(start + dim(arrs[[i]])[3]-1)] <- arrs[[i]]
    start <- start + dim(arrs[[i]])[3]
  }
  b <- as.bugs.array(arr, n.burnin = m$mcmc$burnin, n.thin = m$mcmc$thin)
  return(b)
}





