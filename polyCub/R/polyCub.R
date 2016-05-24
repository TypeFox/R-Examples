################################################################################
### Part of the R package "polyCub".
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Copyright (C) 2009-2013 Sebastian Meyer
### Time-stamp: <[polyCub.R] by SM Sam 06/07/2013 12:52 (CEST)>
################################################################################


#' Wrapper Function for the Various Cubature Methods
#'
#' Instead of calling one of the specific cubature methods of this package, the
#' wrapper function \code{polyCub} may be used together with the \code{method}
#' argument. 
#'
#' @param polyregion a polygonal integration domain.
#' The supported classes depend on the specific method, however, the
#' \code{"\link[spatstat]{owin}"} class from package \pkg{spatstat} works for
#' all methods, as well should a \code{"\link[rgeos:gpc.poly-class]{gpc.poly}"}
#' polygon (but see the comments in \code{help("\link{coerce-methods}")}).
#' @param f two-dimensional function to be integrated.
#' As its first argument the function must take a coordinate matrix, i.e. a
#' numeric matrix with two columns. For the \code{"exact.Gauss"} \code{method},
#' \code{f} is ignored since it is specific to the bivariate normal density.
#' @param method choose one of the implemented cubature methods (partial
#' argument matching is applied), see \code{help("\link{polyCub-package}")}
#' for an overview. Defaults to using the product Gauss cubature
#' implemented in \code{\link{polyCub.SV}}.
#' @param ... arguments of \code{f} or of the specific \code{method}.
#' @param plot logical indicating if an illustrative plot of the numerical
#' integration should be produced.
#' @return The approximated integral of \code{f} over \code{polyregion}.
#' @example inst/examples/polyCub.R
#' @keywords math spatial
#' @family polyCub-methods
#' @export

polyCub <- function (polyregion, f,
                     method = c("SV", "midpoint", "iso", "exact.Gauss"), ...,
                     plot = FALSE)
{
	method <- match.arg(method)
	cl <- match.call()
	cl$method <- NULL
	cl[[1]] <- as.name(paste("polyCub", method, sep="."))
	if (method == "exact.Gauss") cl$f <- NULL
	int <- eval(cl, parent.frame())
	int  #structure(int, method = method)
}
