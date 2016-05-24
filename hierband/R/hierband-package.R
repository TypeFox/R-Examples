#' Convex banding of the covariance matrix using 
#'
#' \code{hierband} is the R package implementing the convex banding approach to covariance estimation of Bien, Bunea, & Xiao (see full reference below).
#' 
#' The package is designed for situations in which there is a large number of variables that have a known ordering and in which it is believed that variables far apart in this ordering have little to no dependence.
#'
#' It is called \code{hierband} (pronounced "hair band") because it makes use of a hierarchical group lasso penalty and provides a banded estimate of the covariance matrix.
#' 
#' Its main functions are \code{\link{hierband}}, \code{\link{hierband.path}}, \code{\link{hierband.cv}}.
#' 
#' The development of this package was supported by National Science Foundation grant DMS-1405746.
#'
#' @references Bien, J., Bunea, F. Xiao, L. (2014) "Convex banding of the covariance matrix." Accepted for publication in JASA.
#' @name hierband-package
#' @author Jacob Bien \email{jbien@@cornell.edu}, Florentina Bunea, Luo Xiao
#' @docType package
NULL
