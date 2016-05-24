#' Screen voxels for a voxelwise smoothing object
#' 
#' Inputs a voxelwise smoothing object as produced by \code{\link{semipar4d}},
#' and outputs an object containing the results for a subset of the voxels.
#' 
#' 
#' @param semi.obj an object of class \code{\link{semipar.mp}}.
#' @param arr4d the 4-dimensional array used to generate the object.
#' @param include a logical matrix indicating which points (or voxels) should
#' be included.
#' @return a modified version of \code{semipar.obj}, with pointwise
#' coefficients (\code{coef} component), pointwise degrees of freedom
#' (\code{pwdf}), pointwise log smoothing parameter (\code{pwlsp}), and
#' pointwise variance estimate (\code{sigma2}) for the points specified by
#' \code{include} only.
#' @author Lei Huang \email{huangracer@@gmail.com} and Philip Reiss
#' \email{phil.reiss@@nyumc.org}
#' @seealso \code{\link{semipar.mp}}
#' @examples
#' 
#' data(test)
#' d4 = test$d4
#' x = test$x
#' vw.obj = semipar4d(d4, formula = ~sf(x), data = data.frame(x = x), lsp=-5:5)
#' 
#' # Include only the first 600 voxels
#' sv = screen.vox(vw.obj, d4, rep(1:0, c(600,400)))
#' @export
screen.vox <-
function(semi.obj, arr4d, include)    {
    has.data = attributes(arr4d)$has.data
    data.inds = which(has.data==TRUE, arr.ind=TRUE)
    if (ncol(semi.obj$coef)!=length(include))  stop("Dimensions don't match!")
    semi.obj$coef = semi.obj$coef[ , include]
    semi.obj$pwdf = semi.obj$pwdf[include]
    semi.obj$pwlsp = semi.obj$pwlsp[include]
    semi.obj$sigma2 = semi.obj$sigma2[include]
    semi.obj$include = has.data
    semi.obj$include[has.data & !is.na(has.data)] = include
    return(semi.obj)
}

