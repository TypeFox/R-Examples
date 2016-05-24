##' This function serves to \sQuote{complete} \code{rms:fastbw}, such that
##' it converts one fitted model to another, the restricted model.
##'
##' @title fastback
##' @param fit A model fit of class \code{rms}
##' @param data The data against which \code{fit} was estimated
##' @param ... Other parameters to be passed to \code{fastbw}
##' @return A fitted model of the same type as \code{fit}, with
##'         regressors chosen by stepwise backward regression
##' @author David C. Norris
##' @keywords regression
##' @export fastback
fastback <- function(fit, data, ...){
  ## N.B. We suppress nuisance warnings "force probably does not work..." from 'fastbw':
  deleted.vars <- rownames(suppressWarnings(fastbw(fit, ...))$result)
  smaller.formula <- remvar(formula(fit), deleted.vars)
  fit$call[[match("formula", names(fit$call))]] <- smaller.formula
  red <- eval(fit$call)
  ## In case a formula with trivial (just Intercept) RHS was estimated,
  ## red$var will be NULL. To permit this special case to be handled
  ## within the general algorithm below, we replace this NULL value
  ## with a 1x1 zero matrix:
  if(is.null(red$var)){
    red$var <- matrix(0, dimnames=as.list(rep(names(red$coefficients), 2)))
    red$info.matrix <- matrix(0)
  }
  ## To return a restricted model more readily comparable with the full model,
  ## we now zero the coefficients and variances of the original full model,
  ## and fill in those parts which correspond to the (fewer) parts of the
  ## restricted model.
  keep <- names(red$coefficients)
  k <- match(keep, names(fit$coefficients))
  ## Generally speaking, named numeric vectors and matrices of the original fit
  ## should be zeroed, and replaced where possible with corresponding elements
  ## of the reduced model...
  for(elem in names(fit)){
    if(is.numeric(fit[[elem]])){
      if(is.vector(fit[[elem]])){
        if(length(fit[[elem]])==length(red[[elem]]))
          fit[[elem]] <- red[[elem]]
        else if(!is.null(names(fit[[elem]]))){ # Assume named vectors of changed length correspond to coefficients
          fit[[elem]] <- 0 * fit[[elem]]
          fit[[elem]][keep] <- red[[elem]][keep] # NB: An error will occur, if this assumption is violated!
        }
      } else if(is.matrix(fit[[elem]])){
        fit[[elem]] <- 0 * fit[[elem]]
        if(!is.null(dimnames(red[[elem]])))
          fit[[elem]][keep,keep] <- red[[elem]][keep,keep]
        else
          fit[[elem]][k,k] <- red[[elem]] # TODO: Ascertain that this ALWAYS produces correct indexing.
      }
    }
  }
  ## ..but there are a few exceptions that we deal with separately:
  fit$fail <- red$fail
  fit$deviance <- red$deviance
  ## When calculating the Wald statistic, print.lrm encounters a 0/0=NaN error
  ## where parameter estimates and the associated variances are both strictly zero.
  ## Therefore, we set the variances of excluded variables to a vanishingly small,
  ## but finite value:
  fit$var[fit$var==0] <- .Machine$double.eps
  ## Return the reduced fit
  eval(fit) # TODO: Find out whether the eval(.) really makes the difference!
}
