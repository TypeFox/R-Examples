#' Voxelwise restricted likelihood ratio tests
#' 
#' A wrapper function for \code{\link{rlrt.mp}} to handle 3D image responses.
#' 
#' 
#' @param arr4d a 4-dimensional response array, where the first 3 dimensions
#' refer to spatial coordinates and the last dimension corresponds to different
#' images.
#' @param x,nbasis,norder,nulldim,loginvsp,get.df,B,P see
#' \code{\link{rlrt.mp}}.
#' @return A massively parallel RLRT object, as produced by
#' \code{\link{rlrt.mp}}.
#' @author Lei Huang \email{huangracer@@gmail.com} and Philip Reiss
#' \email{phil.reiss@@nyumc.org}
#' @seealso \code{\link{plot.rlrt4d}}, \code{\link{rlrt.mp}}
#' @examples
#' 
#' data(test)
#' d4 = test$d4
#' x = test$x
#' rlrtobj = rlrt4d(d4, x, loginvsp = -5:5)
#' plot(rlrtobj, d4, slice=5)

#' \dontrun{rlrtpanel(rlrtobj, d4, x)}
#' @export
rlrt4d <-
function(arr4d, x=NULL, nbasis = 15, norder=4, nulldim=NULL, loginvsp, get.df=FALSE, B=NULL, P=NULL) {
    ndim = length(dim(arr4d))
    has.data <- attributes(arr4d)$has.data
    N = NROW(x)
    dim(arr4d) = c(prod(dim(arr4d)[1:(ndim-1)]), N)
    Y.d = t(arr4d[(as.vector(has.data)), ])
    rm(arr4d)  # i.e., remove from current environment
    rlobj <- rlrt.mp(Y.d, x=x, nbasis = nbasis, norder=norder, nulldim=nulldim, loginvsp=loginvsp, get.df=get.df, B=B, P=P)
    rlobj$call = match.call()
    class(rlobj) = "rlrt4d"
    return(rlobj)
}
