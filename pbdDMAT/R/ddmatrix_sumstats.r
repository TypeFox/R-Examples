#' Basic Summary Statistics
#' 
#' Get basic summary statistics.
#' 
#' The return is on process 0 only.
#' 
#' @param x 
#' numeric distributed matrix
#' @param ...
#' Additional arguments.
#' @param na.rm
#' Handling of NA's.
#' 
#' @return 
#' A single value, owned by all ranks in the MPI communicator.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:16, ncol=4)
#' dx <- as.ddmatrix(x) 
#' 
#' summary(dx)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name ddmatrix-sumstats
#' @rdname ddmatrix-sumstats
NULL



#' @rdname ddmatrix-sumstats
#' @export
setMethod("sum", signature(x="ddmatrix"),
  function(x, ..., na.rm=FALSE)
  {
    # no need to correct for local storage issues
    other <- list(...)
    if (length(other) > 0)
      other <- sum(
        sapply(other, 
          function(i) {
            if (is.ddmatrix(i)) 
              sum(i@Data, na.rm=na.rm) 
            else {
              if (comm.rank()==0)
                sum(i, na.rm=na.rm)
              else
                0
            }
          }
        ), 
      na.rm=na.rm)
    else
      other <- 0
    
    local <- sum(x@Data, na.rm=na.rm) + other
    pbdMPI::allreduce(local, op="sum")
  }
)

# mean, with large chunks taken from base:::mean.default
#' @rdname ddmatrix-sumstats
#' @export
setMethod("mean", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (na.rm) 
        x@Data <- matrix(x@Data[!is.na(x@Data)])
#    if (!is.numeric(trim) || length(trim) != 1L) 
#      comm.stop("'trim' must be numeric of length one")
    if (!base.ownany(x@dim, x@bldim, x@ICTXT))
      n <- 0
    else
    n <- length(x@Data)
    n <- pbdMPI::allreduce(n, op='sum')
#    if (trim > 0 && n) {
#        if (is.complex(x)) 
#            comm.stop("trimmed means are not defined for complex data")
#        if (any(is.na(x))) 
#            return(NA_real_)
#        if (trim >= 0.5) 
#            return(median(x, na.rm = FALSE))
##        lo <- floor(n * trim) + 1
##        hi <- n + 1 - lo
##        x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
#    }
    
    sum(x, na.rm=na.rm) / n
  }
)

# prod
#' @rdname ddmatrix-sumstats
#' @export
setMethod("prod", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      prod <- prod(x@Data, na.rm=na.rm)
    else
      prod <- 1
    pbdMPI::allreduce(prod, op="prod")
  }
)

# min/max
#' @rdname ddmatrix-sumstats
#' @export
setMethod("min", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      min <- min(x@Data, na.rm=na.rm)
    else
      min <- Inf
    pbdMPI::allreduce(min, op="min")
  }
)

#' @rdname ddmatrix-sumstats
#' @export
setMethod("max", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (base.ownany(dim=x@dim, bldim=x@bldim, ICTXT=x@ICTXT))
      max <- max(x@Data, na.rm=na.rm)
    else
      max <- -Inf
    pbdMPI::allreduce(max(x@Data, na.rm=na.rm), op="max")
  }
)




# Implementation of http://www.umiacs.umd.edu/research/EXPAR/papers/3494/node18.html#SECTION00051000000000000000
# finds k'th ordered element of the distributed matrix
# desperately needs to be rewritten in C
dmat.rank_k <- function(vec, k, shouldsort=FALSE)
{
  if (shouldsort)
    vec <- sort(vec, na.last=NA)
  
  # FIXME change to numroc call
  if (length(vec)==1)
    if (vec==0)
      vec <- NA
  
    mdmd <- median(unlist(pbdMPI::allgather(median(vec, na.rm=TRUE))), na.rm=TRUE)

    below <- vec[which(vec <= mdmd)]
    lbelow <- length(below)
    test <- pbdMPI::allreduce(lbelow, op='sum')

    if (test < k){
      vec <- vec[which(vec > mdmd)]
      k <- k - test
      mdmd <- dmat.rank_k(vec=vec, k=k)
    }
    else if (test > k){
      vec <- below
      mdmd <- dmat.rank_k(vec=vec, k=k)
    } else {
    
      if (lbelow==0)
        below <- -Inf
      else if (is.na(below)[1])
        below <- -Inf
      mxbl <- max(below)
      closest <- diff(c(mxbl, mdmd))
      allclosest <- pbdMPI::allreduce(closest, op='min')
      if (allclosest==closest)
        mdmd <- mxbl
      else
        mdmd <- 0
      mdmd <- pbdMPI::allreduce(mdmd, op='sum')
    }
    
  return(mdmd)
}

rank_k <- dmat.rank_k



#' @rdname ddmatrix-sumstats
#' @export
setMethod("median", signature(x="ddmatrix"),
  function(x, na.rm=FALSE)
  {
    if (!na.rm){
      test <- any(is.na(x@Data))
      test <- pbdMPI::allreduce(test, op='max')
      if (test>0)
        return(NA)
    } else
      x@Data <- matrix(x@Data[!is.na(x@Data)])
    lenloc <- length(x@Data)
    if (!base.ownany(x@dim, x@bldim, x@ICTXT))
      lenloc <- 0
    n <- pbdMPI::allreduce(lenloc, op='sum')
    if (n%%2==1)
      ret <- dmat.rank_k(vec=x@Data, k=ceiling(n/2), shouldsort=TRUE)
    else {
      ret1 <- dmat.rank_k(vec=x@Data, k=ceiling(n/2), shouldsort=TRUE)
      ret2 <- dmat.rank_k(vec=x@Data, k=ceiling((n+1)/2), shouldsort=TRUE)
      ret <- mean(c(ret1, ret2))
    }
    
    return( ret )
  }
)

