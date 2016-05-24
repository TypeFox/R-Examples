#' Log-likelihood computation for distance sampling data
#'
#' For a specific set of parameter values, it computes and returns the negative
#' log-likelihood for the distance sampling likelihood for distances that are
#' unbinned, binned and a mixture of both.  The function \code{flnl} is the
#' function minimized using \code{\link{optim}} from within
#' \code{\link{ddf.ds}}.
#'
#' Most of the computation is in \code{flpt.lnl} in which the negative
#' log-likelihood is computed for each observation. \code{flnl} is a wrapper
#' that optionally outputs intermediate results and sums the individual
#' log-likelihood values.
#'
#' \code{flnl} is the main routine that manipulates the parameters using
#' \code{\link{getpar}} to handle fitting of key, adjustment or all of the
#' parameters.  It then calls \code{flpt.lnl} to do the actual computation of
#' the likelihood.  The probability density function for point counts is
#' \code{fr} and for line transects is \code{fx}.
#' \code{fx}=g(x)/mu (where g(x) is the detection function); whereas,
#' f(r)=r*g(r)/mu where mu in both cases is the normalizing constant.  Both
#' functions are in source code file for \code{link{detfct}} and are called from
#' \code{distpdf} and the integral calculations are made with
#' \code{\link{integratepdf}}.
#'
#' @aliases flnl flpt.lnl
#' @param fpar parameter values for detection function at which log-likelihood
#'   should be evaluated
#' @param ddfobj distance sampling object
#' @param misc.options width-transect width (W); int.range-integration range
#'   for observations; showit- 0 to 3 controls level of iteration output; 
#'   integral.numeric-if TRUE integral is computed numerically rather
#'   than analytically
#' @param fitting "key" if only fitting key fct parameters, "adjust" if fitting
#'   adjustment function parameters or "all" to fit both
#' @return negative log-likelihood value at the parameter values specified in
#'   \code{fpar}
#' @note These are internal functions used by \code{\link{ddf.ds}} to fit
#'   distance sampling detection functions.  It is not intended for the user to
#'   invoke these functions but they are documented here for completeness.
#' @author Jeff Laake, David L Miller
#' @seealso \code{\link{flt.var}}, \code{\link{detfct}}
#' @keywords utility
flnl <- function(fpar, ddfobj, misc.options, fitting="all"){

  # During the optimisation we want to make sure that we are keeping the
  # right things constant, so lets do that...
  set.na.pars <- function(par.name,ddfobj,fpar){
    # if the parameters exist
    if(!is.null(ddfobj[[par.name]])){
      # set those we don't want to optimise as NA
      save.pars <- ddfobj[[par.name]]$parameters
      ddfobj[[par.name]]$parameters <- rep(NA,
                                          length(ddfobj[[par.name]]$parameters))
      pars <- getpar(ddfobj)
      fpar[which(is.na(pars))] <- save.pars
    }
    return(list(fpar=fpar,ddfobj=ddfobj))
  }

  if(fitting=="key"){
    setna <- set.na.pars("adjustment",ddfobj,fpar)
    ddfobj <- setna$ddfobj
    fpar <- setna$fpar
  }else if(fitting=="adjust"){
    setna <- set.na.pars("shape",ddfobj,fpar)
    ddfobj <- setna$ddfobj
    fpar <- setna$fpar

    setna <- set.na.pars("scale",ddfobj,fpar)
    ddfobj <- setna$ddfobj
    fpar <- setna$fpar
  }

  #  compute total negative log-likelihood
  lnl <- sum(flpt.lnl(fpar, ddfobj, misc.options))

  # If iteration results are printed, output
  # log-likelihood and parameter values
  if(misc.options$showit >= 3){
    cat("par = ", fpar,"\n")
    cat("lt lnl = ", lnl,   "\n")
  }
  return(lnl)
}
