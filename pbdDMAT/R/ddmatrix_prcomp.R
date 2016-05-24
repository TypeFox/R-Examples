# Large parts of these R wrappers taken from core R prcomp.default and 
# scale.default functions

#' Principal Components Analysis
#' 
#' Performs the principal components analysis.
#' 
#' \code{prcomp()} performs the principal components analysis on the data
#' matrix by taking the SVD. Sometimes core R and pbdDMAT will disagree
#' slightly in what the rotated variables are because of how the SVD is
#' caluclated.
#' 
#' @param x 
#' numeric distributed matrix.
#' @param center 
#' logical value, determines whether or not columns are zero
#' centered
#' @param scale. 
#' logical value, determines whether or not columns are rescaled
#' to unit variance
#' @param retx 
#' logical, indicates whether the rotated variables should be
#' returned
#' @param tol 
#' a value indicating the magnitude below which components should be
#' omitted. (Components are omitted if their standard deviations are less than
#' or equal to \code{tol} times the standard deviation of the first component.)
#' With the default null setting, no components are omitted.  Other settings
#' for tol could be \code{tol = 0} or \code{tol = sqrt(.Machine$double.eps)},
#' which would omit essentially constant components
#' @param ...
#' Ignored.
#' 
#' @return 
#' Returns a list.
#' 
#' @keywords Methods
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' comm.set.seed(diff=T)
#' 
#' x <- ddmatrix("rnorm", 10, 10)
#' 
#' y <- prcomp(x)
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
#' @name ddmatrix-prcomp
#' @rdname ddmatrix-prcomp
NULL



#' @rdname ddmatrix-prcomp
#' @export
setGeneric(name = "prcomp", useAsDefault = stats::prcomp, package="pbdDMAT")



#' @rdname ddmatrix-prcomp
#' @export
setMethod("prcomp", signature(x="ddmatrix"),
  function(x, retx=TRUE, center=TRUE, scale.=FALSE, tol=NULL, ...) 
  {
      x <- scale(x, center = center, scale = scale.)
      cen <- attr(x@Data, "scaled:center")
      sc <- attr(x@Data , "scaled:scale")
      if (any(sc == 0))
          comm.stop("cannot rescale a constant/zero column to unit variance")
      
      s <- svd(x, nu=0)
      s$d <- s$d/sqrt(max(1, nrow(x) - 1))
      if (!is.null(tol)) {
          rank <- max(sum(s$d > (s$d[1L] * tol)), 1)
          if (rank < ncol(x)) {
              s$v <- s$v[, 1L:rank]
              s$d <- s$d[1L:rank]
          }
      }
      r <- list(sdev = s$d, rotation = s$v, 
                center = if (is.null(cen)) FALSE else cen, 
                scale = if (is.null(sc)) FALSE else sc)
      if (retx) r$x <- x %*% s$v
      class(r) <- "prcomp"
      
      return(r)
  }
)


