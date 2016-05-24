#' Extract curve estimates to be clustered 
#'
#' Given a massively parallel smoothing object created by \code{\link{semipar.mp}}, this function extracts an object of
#' class \code{"\link{fd}"} representing the curves that one wishes to cluster using \code{\link{funkmeans}}.
#'
#' @param obj object created by \code{\link{semipar.mp}}.
#' @param term which smooth term to extract (useful if the fitted model includes more than one term).
#' @param intercept logical; if \code{TRUE}, intercept will be added to all coefficients. For simple nonparametric regression
#' this should be done to recover the fitted values. 
#' @return an object of class \code{"\link{fd}"} representing the fitted curves, which can be clustered by 
#' \code{\link{funkmeans}}.
#' @author Ruixin Tan
#' @export
#' @import fda
#' @seealso \code{\link{semipar.mp}}, \code{\link{funkmeans}}
#' @examples
#' # see example for plot.funkmeans

extract.fd = function(obj, term = 1, intercept = (term == 1)) {

    # give error if 'term' doesn't refer to one of the smooth terms
    if(!(term %in% obj$where.sf) )  stop("'term' must refer to one of the smooth terms of 'obj'")

    # use 'term' to extract coefficient matrix 
    temp = obj$list.all[[term]]
    coef = obj$coef[temp$start:temp$end, ]
    if (intercept){
        int = obj$coef[1, ]
        add.int = function(x) x[-1] + x[1]
        coef = apply(rbind(int, coef), 2, add.int)
    }

    # generate fdobj
    fdobj = fd(coef = coef, basisobj = temp$basis)

}
