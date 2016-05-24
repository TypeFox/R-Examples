#' @include generics.R
{}



#' Archetypes object constructor
#'
#' @param object The archetypes; a \eqn{p \times m} matrix, see
#'   \code{\link{parameters}}.
#' @param k The number of archetypes;
#' @param alphas The coefficients; a \eqn{n \times p} matrix, see
#'   \code{\link{coef}}.
#' @param rss The residual sum of squares; see \code{\link{rss.archetypes}}.
#' @param iters The number of iterations to the convergence.
#' @param call The call of the \code{\link{archetypes}} function.
#' @param history If \code{saveHistory} set then an environment with the
#'   archetypes object for each execution step;
#' @param kappas The kappas for each system of linear equations.
#' @param betas The data coefficients; a \eqn{p \times n} matrix.
#' @param zas The temporary archetypes.
#' @param family The archetypes family.
#' @param familyArgs Additional arguments for family blocks.
#' @param residuals The residuals.
#' @param weights The data weights.
#' @param reweights The data reweights.
#' @param scaling The scaling parameters of the data.
#'
#' @return A list with an element for each parameter and class attribute
#'   \code{archetypes}.
#'
#' @family archetypes
#'
#' @export
as.archetypes <- function(object, k, alphas, rss, iters = NULL, call = NULL,
                          history = NULL, kappas = NULL, betas = NULL, zas = NULL,
                          family = NULL, familyArgs = NULL, residuals = NULL,
                          weights = NULL, reweights = NULL, scaling = NULL) {

  return(structure(list(archetypes = object,
                        k = k,
                        alphas = alphas,
                        rss = rss,
                        iters = iters,
                        kappas = kappas,
                        betas = betas,
                        zas = zas,
                        call = call,
                        history = history,
                        family = family,
                        familyArgs = familyArgs,
                        residuals = residuals,
                        weights = weights,
                        reweights = reweights,
                        scaling = scaling),
                   class = c(family$class, 'archetypes')))
}



setOldClass(c("archetypes"))


#' @S3method print archetypes
print.archetypes <- function(x, full = TRUE, ...) {
  if ( full ) {
    cat('Archetypes object\n\n')
    cat(paste(deparse(x$call), collapse = '\n'), '\n\n')
  }

  cat('Convergence after', x$iters, 'iterations\n')
  cat('with RSS = ', rss(x), '.\n', sep = '')
}



#' Return fitted data
#'
#' Returns the approximated data.
#'
#' @param object An \code{archetypes} object.
#' @param ... Ignored.
#' @return Matrix with approximated data.
#' @method fitted archetypes
#' @rdname fitted
#'
#' @importFrom stats fitted
#' @S3method fitted archetypes
fitted.archetypes <- function(object, ...) {
  t(t(object$archetypes) %*% t(object$alphas))
}



#' Return fitted archetypes
#'
#' @param object An \code{archetypes} object.
#' @param ... Ignored.
#' @return Matrix with \eqn{k} archetypes.
#'
#' @aliases parameters-methods
#' @aliases parameters,archetypes-method
#'
#' @import methods
#' @importFrom modeltools parameters
#' @exportMethod parameters
#' @rdname parameters
setMethod('parameters', signature = c(object = 'archetypes'),
function(object, ...) { 
  object$archetypes
})



#' Return coefficients
#'
#' @param object An \code{archetypes} object.
#' @param type Return alpha or beta coefficients.
#' @param ... Ignored.
#' @return Coefficient matrix.
#' @method coef archetypes
#' @rdname coef
#'
#' @importFrom stats coef
#' @S3method coef archetypes
coef.archetypes <- function(object, type = c('alphas', 'betas'), ...) {
  type <- match.arg(type)
  object[[type]]
}



#' Return residuals
#'
#' @param object An \code{archetypes} object.
#' @param ... Ignored.
#' @return Matrix with residuals.
#' @method residuals archetypes
#' @rdname residuals
#'
#' @importFrom stats residuals
#' @S3method residuals archetypes
residuals.archetypes <- function(object, ...) {
  object$residuals
}



#' Return residual sum of squares
#'
#' @param object An \code{archetypes} object.
#' @param type Return scaled, single or global RSS.
#' @param ... Ignored.
#' @return Residual sum of squares.
#' @method rss archetypes
#' @rdname rss
#'
#' @S3method rss archetypes
rss.archetypes <- function(object, type = c('scaled', 'single', 'global'), ...) {
  type <- match.arg(type)
  resid <- residuals(object)

  switch(type,
         scaled = object$rss,
         single = apply(resid, 1, object$family$normfn),
         global = object$family$normfn(resid) / nrow(resid))
}



#' Return weights
#'
#' @param object An \code{archetypes} object.
#' @param type Return global weights (weighted archetypes) or
#'   weights calculated during the iterations (robust archetypes).
#' @param ... Ignored.
#' @return Vector of weights.
#' @method weights archetypes
#' @rdname weights
#'
#' @importFrom stats weights
#' @S3method weights archetypes
weights.archetypes <- function(object, type = c('weights', 'reweights'), ...) {
  type <- match.arg(type)
  object[[type]]
}



#' Return kappa
#'
#' @param z An \code{archetypes} object.
#' @param ... Ignored.
#' @return A vector of kappas.
#' @rdname kappa
#'
#' @method kappa archetypes
#' @S3method kappa archetypes
kappa.archetypes <- function(z, ...) {
  return(z$kappas)
}



#' Return number of archetypes
#'
#' @param object An \code{archetypes} object.
#' @param ... Ignored.
#' @return Number of archetypes.
#' @rdname nparameters
#'
#' @method nparameters archetypes
#' @S3method nparameters archetypes
nparameters.archetypes <- function(object, ...) {
  return(object$k)
}



#' Predict method for archetypal analysis fits
#'
#' This method produces predicted alpha coefficients for new data.
#'
#' @param object An \code{archetypes} object; currently only
#'   \code{\link[=archetypesFamily]{original}}-family objects.
#' @param newdata A data frame with data for which to
#'   predict the alpha coefficients.
#' @param ... Ignored.
#' @return The predict alpha coefficients.
#' @rdname predict
#'
#' @method predict archetypes
#' @S3method predict archetypes
predict.archetypes <- function(object, newdata, ...) {
  stopifnot(object$family$which == "original")

  scale <- object$scaling

  ## HACK: use blocks!
  x <- t(newdata)
  x <- x - scale$mean
  x <- x / scale$sd
  x <- object$family$dummyfn(x, ...)

  zs <- t(parameters(object))
  zs <- zs - scale$mean
  zs <- zs / scale$sd
  zs <- rbind(zs, 200)

  alphas <- matrix(NA, ncol = ncol(x), nrow = ncol(coef(object)))
  alphas <- object$family$alphasfn(alphas, zs, x)

  t(alphas)
}
