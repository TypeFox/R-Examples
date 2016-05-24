#' Covariance and Correlation
#' 
#' \code{cov()}
#' and \code{var()} form the variance-covariance matrix.  \code{cor()} forms
#' the correlation matrix.  \code{cov2cor()} scales a covariance matrix into a
#' correlation matrix.
#' 
#' \code{cov()} forms the variance-covariance matrix. Only
#' \code{method="pearson"} is implemented at this time.
#' 
#' \code{var()} is a shallow wrapper for \code{cov()} in the case of a
#' distributed matrix.
#' 
#' \code{cov2cor()} scales a covariance matrix into a correlation matrix.
#' 
#' @param x,y,V 
#' numeric distributed matrices.
#' @param na.rm 
#' logical, determines whether or not \code{NA}'s should be dealth
#' with.
#' @param use 
#' character indicating how missing values should be treated.
#' Acceptable values are the same as \code{R}'s, namely "everything",
#' "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param method 
#' character argument indicating which method should be used to
#' calculate covariances. Currently only "spearman" is available for
#' \code{ddmatrix}.
#' 
#' @return 
#' Returns a distributed matrix.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' x <- ddmatrix("rnorm", nrow=3, ncol=3), bldim=2
#' 
#' cv <- cov(x)
#' print(cv)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name covariance
#' @rdname covariance
NULL



#' @rdname covariance
#' @export
setMethod("cov", signature(x="ddmatrix"),
function (x, y = NULL, use = "everything", method = "pearson") 
  {
    yexists <- !is.null(y)
    
    if (yexists){
      if (!is.ddmatrix(y))
        comm.stop("Error : 'y' must be a distributed matrix")
      else if (x@dim[1] != y@dim[1])
        comm.stop("Error : incompatible dimensions")
    }
    
    if (use=="all.obs"){
      anyna <- FALSE
      if (any(is.na(x)))
        anyna <- TRUE
      if (yexists)
        if (any(is.na(y)))
          anyna <- TRUE
      
      if (anyna)
        comm.stop("Error : missing observations in cov")
    }
    if (use=="complete.obs"){
      if (yexists){
        narows <- unique(which(is.na(rowSums(x) + rowSums(y))))
        lnarows <- length(narows)
        if (lnarows > 0) {
          if (lnarows < x@dim[1]){
            x <- x[-narows, ]
            y <- y[-narows, ]
          } 
          else 
            comm.stop("Error : no complete element pairs")
        }
      } 
      else
        x <- na.exclude(x)
    } 
    else if (use=="na.or.complete")
      comm.stop("Error : na.or.complete not yet implemented")
    else if (use=="pairwise.complete.obs")
      comm.stop("Error : pairwise.complete.obs not yet implemented")
    else if (use!="everything")
      comm.stop("Error : invalid 'use' argument")
    
    method <- pbdMPI::comm.match.arg(method, c("pearson", "kendall", "spearman"))
    if (method == "pearson") {
#########################################################
      cntr <- dmat.clmn(x, na.rm=FALSE)
      if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
        x@Data <- base::scale(x@Data, center=cntr, scale=FALSE)
#      x <- scale(x, scale=FALSE)
#########################################################
      if (is.null(y))
        ret <- crossprod(x=x) / max(1, nrow(x) - 1)
      else {
        scale(x=y, center=TRUE, scale=FALSE)
        ret <- t(x) %*% y / (nrow(x) - 1)
      }
    }
    else 
      comm.stop("Error : Other methods not yet implemented")
    
    if (x@dim[1] == 1)      
      ret[, ] <- NA
    
    return( ret )
  }
)



# Much of this wrapper taken from core R's var function
#' @rdname covariance
#' @export
setMethod("var", signature(x="ddmatrix"),
function (x, y = NULL, na.rm = FALSE, use) 
  {
    if (missing(use)) {
      if (na.rm) 
        use <- "na.or.complete"
      else 
        use <- "everything"
    }
    
    na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", "everything", "na.or.complete"))
    if (is.na(na.method)) 
        comm.stop("invalid 'use' argument")
    
    ret <- cov(x, y, na.method, FALSE)
    
    return( ret )
  }
)



# experimental
#' @rdname covariance
#' @export
setMethod("cor", signature(x="ddmatrix"),
function (x, y = NULL, use = "everything", method = "pearson") 
  {
    if (method != "pearson")
      comm.stop("Not yet implemented")
    
    if (use != "everything")
      comm.stop("Not yet implemented")
    
    x <- scale(x=x, center=TRUE, scale=TRUE)
    
    if (!is.null(y)){
      yscaled <- scale(x=y, center=TRUE, scale=TRUE)
      ret <- t(x) %*% y / (nrow(x) - 1)
    }
    else {
      ret <- crossprod(x=x) / max(1, nrow(x) - 1)
    }
    
    return( ret )
  }
)



#' @rdname covariance
#' @export
setMethod("cov2cor", signature(V="ddmatrix"),
function(V)
  {
    d <- sqrt(1/diag(V))
    
    r <- V@Data
    descv <- base.descinit(dim=V@dim, bldim=V@bldim, ldim=V@ldim, ICTXT=V@ICTXT)
    
    r <- base.pdsweep(x=r, descx=descv, vec=d, MARGIN=1L, FUN="*")
    r <- base.pdsweep(x=r, descx=descv, vec=d, MARGIN=2L, FUN="*")
    
    V@Data <- r
    
    return( V )
  }
)



