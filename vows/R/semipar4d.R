#' Massively parallel semiparametric regression for 4-dimensional data
#' 
#' This is a wrapper function for \code{\link{semipar.mp}} to handle 3D image
#' responses.
#' 
#' 
#' @param arr4d a 4-dimensional response array, where the first 3 dimensions
#' refer to spatial coordinates and the last dimension corresponds to different
#' images.
#' @param
#' formula,lsp,data,range.basis,knots,rm.constr,random,store.reml,store.fitted
#' see \code{\link{semipar.mp}}.
#' @return An object of class \code{"\link{semipar.mp}"}, with two changes. (1)
#' If \code{store.fitted = TRUE}, the fitted values are given as a
#' 4-dimensional array. (2) A \code{call} component is included.
#' @author Yin-Hsiu Chen \email{enjoychen0701@@gmail.com} and Philip Reiss
#' \email{phil.reiss@@nyumc.org}
#' @seealso \code{\link{semipar.mp}}
#' @examples
#' 
#' data(test)
#' d4 = test$d4
#' x = test$x
#' semi.obj = semipar4d(d4, ~sf(x), lsp=-5:5, data=data.frame(x = x))
#' plot(semi.obj, which.vox = 4)
#' @export
semipar4d <-
function(arr4d, formula, lsp, data, range.basis = NULL, knots = "quantile", rm.constr = FALSE, random = NULL, store.reml = FALSE, store.fitted = FALSE)   {
    dim.4d = dim(arr4d)
    has.data <- attributes(arr4d)$has.data
    N = NROW(data)
    dim(arr4d) = c(prod(dim(arr4d)[1:3]), N)
    Y.d = t(arr4d[(as.vector(has.data)), ])
    rm(arr4d) # i.e., remove from current environment  
    
    Y.fit = semipar.mp(formula = formula, Y = Y.d, lsp = lsp, data = data, range.basis = range.basis, knots = knots, rm.constr = rm.constr, store.reml = store.reml, store.fitted = store.fitted)
    if (store.fitted) {
        fit.value = array(NA, dim=dim.4d)
        for (j in 1:N) fit.value[ , , ,j][has.data] = Y.fit$fitted[j,]
        Y.fit$fitted = fit.value
    }
    Y.fit$call = match.call()
    return(Y.fit)
}
 
