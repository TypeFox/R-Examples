#' Functional k-means clustering for parallel smooths
#' 
#' This function performs k-means clustering for curve estimates corresponding
#' to each of a 3D grid of points. For example, when scatterplot smoothing is
#' performed at each of a grid of brain voxels as in Reiss et al. (2014), this
#' function can be used to cluster the obtained smooths.
#' 
#' The functional clustering algorithm consists of performing (i) functional
#' principal component analysis of the curve estimates or their derivatives,
#' followed by (ii) k-means clustering of the functional PC scores (Tarpey and
#' Kinateder, 2003).
#' 
#' @param fdobj a functional data object, of class \code{"\link[fda]{fd}"},
#' defining the set of curves being clustered.
#' @param deriv which derivative of the curves should be clustered.  If
#' \code{0}, the curves themselves are clustered; if \code{1} (the default),
#' their first derivatives are clustered, a natural way to assign curves of
#' similar shape to the same cluster.
#' @param lambda smoothing parameter for functional PCA as implemented by
#' \code{\link[fda]{pca.fd}}.
#' @param ncomp number of functional principal components.
#' @param centers number of clusters.
#' @param nstart number of randomly chosen sets of initial centers used by the
#' \code{\link[stats]{kmeans}} function.
#' @param store.fdobj logical: Should the input fd object be stored in the
#' output? May wish to set to FALSE for large sets of smooths.
#' @return An object of class "funkmeans", which is a list with elements:
#' \item{cluster, centers, withinss, tots, tot.withinss, betweenness, size}{see
#' \code{\link[stats]{kmeans}}.} \item{basis,coef}{basis object and coefficient
#' matrix defining the functional data object (see \code{\link[fda]{fd}}) for
#' the curves that are clustered.} \item{fpca}{functional principal components
#' object, output by \code{\link[fda]{pca.fd}}.} \item{R2}{proportion of
#' variance explained by the k clusters.}
#' @author Philip Reiss \email{phil.reiss@@nyumc.org}, Lei Huang
#' \email{huangracer@@gmail.com} and Lan Huo
#' @seealso \code{\link{funkmeans4d}}, \code{\link{funkpanel}}
#' @references Alexander-Bloch, A. F., Reiss, P. T., Rapoport, J., McAdams, H.,
#' Giedd, J. N., Bullmore, E. T., and Gogtay, N. (2014). Abnormal cortical
#' growth in schizophrenia targets normative modules of synchronized
#' development. \emph{Biological Psychiatry}, in press.
#' 
#' Reiss, P. T., Huang, L., Chen, Y.-H., Huo, L., Tarpey, T., and Mennes, M.
#' (2014). Massively parallel nonparametric regression, with an application to
#' developmental brain mapping. \emph{Journal of Computational and Graphical
#' Statistics}, \emph{Journal of Computational and Graphical Statistics},
#' 23(1), 232--248.
#' 
#' Tarpey, T., and Kinateder, K. K. J. (2003).  Clustering functional data.
#' \emph{Journal of Classification}, 20, 93--114.
#' @examples
#' 
#' # See example for funkpanel
#' @export
#' @import fda 
funkmeans = function(fdobj, deriv = 1, lambda = 0, ncomp, centers, nstart = 10, store.fdobj=TRUE) {   
    deriv.fdobj = deriv.fd(fdobj, deriv)
    harmfdPar = fdPar(deriv.fdobj)
    harmfdPar$lambda = lambda
    fpca.obj = pca.fd(deriv.fdobj, nharm = ncomp, harmfdPar)
    km.obj = kmeans(fpca.obj$scores, centers = centers, nstart = nstart, iter.max=100)
    km.obj$basis = fdobj$basis
    km.obj$coef = fdobj$coef
    km.obj$fpca = fpca.obj
    km.obj$R2 = (1 - km.obj$tot.withinss / km.obj$totss) * sum(km.obj$fpca$varprop) 
    if (store.fdobj) km.obj$fdobj = fdobj
    class(km.obj) = "funkmeans"
    km.obj
}
