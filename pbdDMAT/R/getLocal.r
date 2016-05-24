#' getLocal
#' 
#' @param x
#' A distributed matrix.
#' @param gi,gj
#' Global row and column indices, respectively.
#' @param all.rank
#' Logical; if \code{TRUE}, then all processes will hold the desired
#' value on exit.  Otherwise, only the process who owns the local
#' value returns this value, while every other process returns \code{NULL}.
#' @param gridinfo
#' An optional parameter; each local data lookup requires the data
#' contained in \code{gridinfo(ICTXT(x))}.  So you may specify it
#' yourself (and if you are making many function calls, this is
#' preferable performance-wise), or the lookup will be performed
#' for you.
#' 
#' @return 
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT)
#' init.grid()
#' 
#' x <- ddmatrix(1:100, 10, bldim=c(2, 2))
#' 
#' val <- getLocal(x, 5, 1)
#' comm.print(val, all.rank=TRUE)
#' 
#' val <- getLocal(x, 5, 1, all.rank=FALSE)
#' comm.print(val, all.rank=TRUE)
#' 
#' finalize()
#' }
#' 
#' @export
getLocal <- function(x, gi, gj, all.rank=TRUE, gridinfo)
{
  if (class(x) != "ddmatrix")
    comm.stop("Argument 'x' must be of type 'ddmatrix'")
  
  if (missing(gridinfo))
    gridinfo <- base.blacs(ICTXT(x))
  
  g <- g2lcoord(dim(x), bldim(x), gi, gj, gridinfo)
  
  if (length(g) == 2L)
    ret <- x@Data[g[1L], g[2L]]
  
  
  if (all.rank)
  {
    if (length(g) == 1L)
      ret <- 0
    
    ret <- allreduce(ret)
  }
  else if (length(g) == 1)
    ret <- invisible()
  
  return(ret)
}
