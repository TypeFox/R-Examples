#' Defining smooth functions in semiparametric model formulae
#' 
#' This function is called by \code{\link{semipar.mp}} to define B-spline
#' smooths.
#' 
#' 
#' @param argvals a vector or matrix of covariates.
#' @param effect predictor whose effect varies with respect to \code{argvals}.
#' E.g., if the effect of \code{diagnosis} varies with \code{age}, use
#' \code{sf(age, effect = diagnosis)}. Similar to argument \code{by} in
#' \code{\link[mgcv]{s}}.
#' @param k number of B-spline basis functions.
#' @param norder order of B-splines: the default, \code{4}, gives cubic
#' B-splines.
#' @param pen.order order of the penalty, i.e., of the derivative defining the
#' penalty.
#' @param range.basis a numeric vector of length 2 defining the interval over
#' which the B-spline basis is created. If \code{NULL}, set to the range of the
#' variable.
#' @param knots knots placement method for B-spline smoothing. The default,
#' "quantile", places the knots at equally spaced quantiles of the data;
#' "equispaced" gives equally spaced knots.
#' @author Yin-Hsiu Chen \email{enjoychen0701@@gmail.com} and Philip Reiss
#' \email{phil.reiss@@nyumc.org}
#' @export
sf <- function(argvals, effect=NULL, k = 10, norder = 4, pen.order = 2, range.basis = NULL, knots = "quantile") {
    if (is.null(range.basis)) range.basis = range(argvals)
    if (knots == "quantile")	basis = create.bspline.basis(range.basis, breaks = quantile(argvals, seq(0,1,length.out=k-norder+2)), norder=norder)
    else  if (knots == "equispaced") basis = create.bspline.basis(range.basis, norder=norder, nbasis = k)
	modmat = eval.basis(argvals, basis) 
	if (!is.null(effect)) modmat = diag(effect) %*% modmat
    penmat = getbasispenalty(basis, pen.order)   
    constraint = if (is.null(effect)) colSums(modmat) else NULL  
    
    sf.out = list(basis = basis, modmat = modmat, penmat = penmat, 
                  constraint = constraint, argvals = argvals, effect = effect, k = k, 
                  norder = norder, pen.order = pen.order)
    class(sf.out) = "sf"
    sf.out
}

