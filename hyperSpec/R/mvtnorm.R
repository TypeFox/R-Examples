.rmmvnorm <- function (n, mean, sigma) {
  .group <- rep.int (seq_along (n), n)

   ## make indices so that pooled or individual covariance matrices can be used.
  if (length (dim (sigma)) == 3L)
    isigma <- seq_len (dim (sigma) [3])
  else {
    isigma <- rep (1L, nrow (mean))
    dim (sigma) <- c (dim (sigma), 1L)
  }

  tmp <- matrix (NA_real_, sum (n), ncol (mean))
  for (i in seq_along (n))
    tmp [.group == i,] <- mvtnorm::rmvnorm (n [i], mean [i,], sigma [,, isigma [i]])

  attr (tmp, "group") <- .group

  tmp
}

##' @export
##' @name rmmvnorm
setGeneric ("rmmvnorm", .rmmvnorm)


##' Multivariate normal random numbers
##'
##' Interface functions to use \code{\link[mvtnorm]{rmvnorm}} for
##' \code{\link[hyperSpec]{hyperSpec-class}} objects.
##'
##' The \code{mvtnorm} method for hyperSpec objects supports producing multivariate normal data for
##' groups with different mean but common covariance matrix, see the examples.
##'
##' @param n vector giving the numer of cases to generate for each group
##' @param mean matrix with mean cases in rows
##' @param sigma common covariance matrix or array (\code{ncol (mean)} x \code{ncol (mean)} x \code{nrow (mean)}) with individual covariance matrices for the groups.
##' @export
##' @seealso \code{\link[mvtnorm]{rmvnorm}}
##'
##' \code{\link[hyperSpec]{cov}} and \code{\link[hyperSpec]{pooled.cov}} about calculating  covariance of hyperSpec objects.
##' @rdname rmmvnorm
##' @importFrom mvtnorm rmvnorm
##' @aliases rmmvnorm rmmvnorm,hyperSpec-method
##' @docType methods
##' @examples
##' ## multiple groups, common covariance matrix
##' 
##' pcov <- pooled.cov (chondro, chondro$clusters)
##' rnd <- rmmvnorm (rep (10, 3), mean = pcov$mean, sigma = pcov$COV)
##' plot (rnd, col = rnd$.group)

setMethod ("rmmvnorm", signature (n = "numeric", mean = "hyperSpec", sigma = "matrix"),
           function (n, mean, sigma){
             tmp <- .rmmvnorm (n, mean@data$spc, sigma)

           data <- mean [attr (tmp, "group"),, drop = FALSE]
           if (hy.getOption ("gc")) gc ()
           data@data$spc <- tmp
           if (hy.getOption ("gc")) gc ()
           data$.group <- attr (tmp, "group")
           if (hy.getOption ("gc")) gc ()
           data
           })

##' @rdname rmmvnorm
##' @export
setMethod ("rmmvnorm", signature (n = "numeric", mean = "hyperSpec", sigma = "array"),
           function (n, mean, sigma){
             tmp <- .rmmvnorm (n, mean@data$spc, sigma)

           data <- mean [attr (tmp, "group"),, drop = FALSE]
           if (hy.getOption ("gc")) gc ()
           data@data$spc <- tmp
           if (hy.getOption ("gc")) gc ()
           data$.group <- attr (tmp, "group")
           if (hy.getOption ("gc")) gc ()
           data
           })

##' @rdname rmmvnorm
##' @export
setMethod ("rmmvnorm", signature (n = "numeric", mean = "matrix", sigma = "matrix"),
           .rmmvnorm)

##' @rdname rmmvnorm
##' @export
setMethod ("rmmvnorm", signature (n = "numeric", mean = "matrix", sigma = "array"),
           .rmmvnorm)




## produces matrices instead of hyperSpec objects. 
## mapply (rmvnorm, n = 1:3, mean = pcov$mean, MoreArgs= list (sigma = pcov$COV), SIMPLIFY = FALSE))
