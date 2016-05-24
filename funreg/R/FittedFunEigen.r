#' @title fitted method for funeigen object 
#' @description Returns fitted values for a \code{funeigen} object.  
#' @details A \code{funeigen} object represents a principal component analysis
#' of irregular longitudinal data, following the method used by Goldsmith et al. (2011).  
#' @param object A  \code{funeigen} object.
#' @param type A character string, one of the following: \code{functions},
#'  \code{eigenfunctions}, \code{loadings}, \code{eigenvalues}, \code{mean},
#'   \code{centered}, \code{covariance}, \code{noise.variance}, 
#'   \code{midpoints}. These are the constructs for which fitted values can be returned.
#' @param ... Other optional arguments which may be passed from 
#' other methods but ignored by this one.
#' @return A matrix or vector containing the appropriate fitted values.  What is 
#' returned depends on the \code{type} parameter.  \code{functions} gives the fitted 
#' values of the smooth latent x(t) functions at a grid of time points.  
#' \code{eigenfunctions} gives the estimated eigenfunctions at each time point.  
#' \code{loadings} gives the loading of each subject on each estimated eigenfunction. 
#'  \code{mean} gives the mean value for the smooth latent x(t) functions.  
#'  \code{centered} gives the centered x(t) functions (the estimated function 
#'  subtracting the mean function) .  \code{covariance} gives the estimated 
#'  covariance matrix of x(s) and x(t) on a grid of time points s and t.
#'  \code{noise.variance} gives the estimated measurement error variance on 
#'  the x(t) functions.  \code{midpoints} gives the time points for the grid, on
#'  which \code{functions}, \code{mean}, \code{centered}, and \code{covariance}
#'  are defined; they are viewed as midpoints of bins of observation times (see 
#'  Goldsmith et al., 2011). 
#' @references  Goldsmith, J., Bobb, J., Crainiceanu, C. M., Caffo, B., and Reich, D.
#'    (2011). Penalized functional regression. Journal of Computational
#'    and Graphical Statistics, 20(4), 830-851. 
#' @S3method fitted funeigen  
#' @method fitted funeigen  
#' @export
fitted.funeigen <- function(object,type="functions",...) {
    type <- tolower(type);
    answer <- NULL;
    if (type=="functions") {
        return(object$x.hat.by.bin);
    }
    if (type=="eigenfunctions") {
        return(object$psi);
    }
    if (type=="loadings") {
        return(object$C);
    }
    if (type=="eigenvalues") {
        return(object$lambda);
    }
    if (type=="mean") {
        return(object$mu.x.by.bin);
    }
    if (type=="centered") {
        return(object$centered.x.hat.by.bin);
    }
    if (type=="covariance") {
        return(object$smoothed.cov.mat.between.bins);
    }
    if (type=="noise.variance") {
        return(object$estimated.noise.variance);
    }
    if (type=="midpoints") {
        return(object$bin.midpoints);
    }
    return(NULL);
}