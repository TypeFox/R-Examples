#' @name get.LML
#' @aliases get.LML
#' @title Log marginal likelihood of model
#' @description Estimate log marginal likelihood of model
#' @param counts vector of counts of the number of proposals that were
#' generated before accepting a draw.  Length of vector is equal to
#' the number of draws from the posterior.  If the first proposal for
#' a particular posterior draw is accepted, that count is a 1.
#' @param log.phi Numeric vector of draws of log.phi from the proposal draws.
#' @param post.mode The posterior mode.
#' @param fn.dens.post Function that returns the log posterior
#' density.  Function should take the parameter vector as the first
#' argument.  Additional arguments are passed as ...
#' @param fn.dens.prop Function that returns the log density of the
#' proposal distribution. The first argument of the function should
#' take either a vector	or a matrix.  If the argument is a matrix,
#' each row is considered a sample.  Additional parameters are passed
#' as a list, prop.params.
#' @param prop.params Object (list or vector) to be passed to both
#' fn.dens.prop and fn.draw.prop.Contains parameters for the proposal
#' distribution.  See details.
#' @param ... Additional parameters to be passed to fn.dens.post
#' @return The estimate log marginal likelihood of the model.
#' @export
get.LML <- function(counts, log.phi, post.mode, fn.dens.post,
                    fn.dens.prop, prop.params, ...) {

  if ( any(!is.finite(counts)) | any(counts<=0) ) {
    stop ("Error in get.LML:  all counts must be finite and positive")
  }

  if ( any(!is.finite(log.phi)) | any(log.phi>0) ) {
    stop ("Error in get.LML:  all values of log.phi must be non-positive")
  }

  log.c1 <- fn.dens.post(post.mode, ...)
  log.c2 <- fn.dens.prop(post.mode, prop.params)
  ar <- 1/mean(counts)
  M <- length(log.phi)
  ord <- order(log.phi,decreasing=TRUE)
  v <- log.phi[ord]
  max.v <- max(v)
  z <- sum((2*seq(1,M)-1) * exp(v-max.v))

  LL <- max.v - log.c2 + log.c1 - log(ar) - 2*log(M) + log(z)
  return(as.numeric(LL))
}

