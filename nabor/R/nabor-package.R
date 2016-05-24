#' Wrapper for libnabo K Nearest Neighbours C++ library
#' 
#' R package \bold{nabor} wraps the 
#' \href{https://github.com/ethz-asl/libnabo}{libnabo} library, a fast K Nearest
#' Neighbour library for low-dimensional spaces written in templated C++. The 
#' package provides both a standalone function (see \code{\link{knn}} for basic 
#' queries along an option to produce an object containing the k-d tree search
#' (see \code{\link{WKNN}}) structure when making multiple queries against the
#' same target points.
#' 
#' libnabo uses the same approach as the ANN library (wrapped in R package 
#' \code{RANN}) but is generally faster and with a smaller memory footprint. 
#' Furthermore since it is templated on the underlying scalar type for 
#' coordinates (among other things), we have provided both float and double 
#' coordinate implementations of the classes wrapping the search tree 
#' structures. See the github repository and Elsenberg et al paper below for 
#' details.
#' @name nabor-package
#' @aliases nabor
#' @useDynLib nabor
#' @import Rcpp methods
#' @references Elseberg J, Magnenat S, Siegwart R and Nuechter A (2012). 
#'   "Comparison of nearest-neighbor-search strategies and implementations for 
#'   efficient shape registration." _Journal of Software Engineering for 
#'   Robotics (JOSER)_, *3*(1), pp. 2-12. ISSN 2035-3928.
#' @seealso \code{\link{knn}}, \code{\link{WKNN}}
NULL

#' List of 3 matrices containing 3D points from Drosophila neurons
#'
#' This R list contains 3 skeletonized \emph{Drosophila} Kenyon cells as
#' \code{dotprops} objects. Original data is due to Chiang et al. 2011, who have
#' generously shared their raw data at \url{http://flycircuit.tw}. Image
#' registration and further processing was carried out by Greg Jefferis.
#' @name kcpoints
#' @docType data
#' @references [1] Chiang A.S., Lin C.Y., Chuang C.C., Chang H.M., Hsieh C.H.,
#'   Yeh C.W., Shih C.T., Wu J.J., Wang G.T., Chen Y.C., Wu C.C., Chen G.Y.,
#'   Ching Y.T., Lee P.C., Lin C.Y., Lin H.H., Wu C.C., Hsu H.W., Huang Y.A.,
#'   Chen J.Y., et al. (2011). Three-dimensional reconstruction of brain-wide
#'   wiring networks in Drosophila at single-cell resolution. Curr Biol 21 (1),
#'   1--11.
NULL
