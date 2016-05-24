
#' Memory-efficient k-means cluster analysis
#'
#' @description k-means cluster analysis without the memory overhead, and 
#' possibly in parallel using shared memory.
#' @param x a \code{\link[bigmemory]{big.matrix}} object.
#' @param centers a scalar denoting the number of clusters, or for k clusters, 
#' a k by \code{ncol(x)} matrix.
#' @param iter.max the maximum number of iterations.
#' @param nstart number of random starts, to be done in parallel if there 
#' is a registered backend (see below).
#' @param dist the distance function. Can be "euclid" or "cosine".
#' @return An object of class \code{kmeans}, just as produced by 
#' \code{\link{kmeans}}.
#' @importFrom foreach foreach %dopar% getDoParName registerDoSEQ 
#' getDoParWorkers
#' @importFrom bigmemory as.big.matrix big.matrix describe
#' @details The real benefit is the lack of memory overhead compared to the 
#' standard \code{\link{kmeans}} function.  Part of the overhead from 
#' \code{kmeans()} stems from the way it looks for unique starting 
#' centers, and could be improved upon.  The \code{bigkmeans()} function 
#' works on either regular \R \code{matrix} objects, or on \code{big.matrix} 
#' objects.  In either case, it requires no extra memory (beyond the data, 
#' other than recording the cluster memberships), whereas \code{kmeans()} 
#' makes at least two extra copies of the data.  And \code{kmeans()} is even 
#' worse if multiple starts (\code{nstart>1}) are used.  If \code{nstart>1} 
#' and you are using \code{bigkmeans()} in parallel, a vector of cluster 
#' memberships will need to be stored for each worker, which could be 
#' memory-intensive for large data.  This isn't a problem if you use are running
#' the multiple starts sequentially.
#'
#' Unless you have a really big data set (where a single run of 
#' \code{\link{kmeans}} not only burns memory but takes more than a few 
#' seconds), use of parallel computing for multiple random starts is unlikely 
#' to be much faster than running iteratively.
#'
#' Only the algorithm by MacQueen is used here.
#' 
#' @note A comment should be made about the excellent package \pkg{foreach}.  By
#' default, it provides \code{\link[foreach]{foreach}}, which is used
#' much like a \code{for} loop, here over the \code{nstart}
# random starting points.  Even so, there are efficiencies, doing a comparison
# of each result to the previous best result (rather than saving everything 
#' and doing a final comparison of all results).
#' 
#' When a parallel backend has been registered (see packages \pkg{doSNOW}, 
#' \pkg{doMC}, and \pkg{doMPI}, for example), \code{bigkmeans()} automatically 
#' distributes the \code{nstart} random starting points across the available 
#' workers.  This is done in shared memory on an SMP, but is distributed on 
#' a cluster *IF* the \code{big.matrix} is file-backed.  If used on a cluster 
#' with an in-RAM \code{big.matrix}, it will fail horribly.  We're considering 
#' an extra option as an alternative to the current behavior.
#' @export
bigkmeans <- function(x, centers, iter.max = 10, nstart = 1, dist='euclid') {

  if (is.null(getDoParName())) {
    registerDoSEQ() # A little hack to avoid the foreach warning 1st time.
  }

  dist_calc = 0
  if (dist=='euclid') {
      dist_calc = 0
  } else if (dist=='cosine') {
      dist_calc = 1
  } else {
      stop("'euclid' or 'cosine' are valid. Check your argument.\n")
  }

  ################################################################
  # This function is used to construct a list of length nstart
  # of centers so that there are no duplicates.  If this dies,
  # the user probably shouldn't be using k-means.
  getcenters <- function(x, k, nstart) {
    n <- nrow(x)
    centers <- list(x[sample(1:n, k),,drop=FALSE])
    nchecks <- 1000
    if (k<=10) nchecks <- 10 + 2^k
    for (ii in 1:nchecks) {
      if (any(duplicated(centers[[length(centers)]]))) {
        centers[[length(centers)]] <- x[sample(1:n, k),,drop=FALSE]
      } else break;
    }
    if (any(duplicated(centers[[length(centers)]]))) {
      stop("Having trouble finding non-duplicated centers.\n")
    }
    if (nstart>1) {
      for (i in 2:nstart) {
        centers[[length(centers)+1]] <- x[sample(1:n, k),,drop=FALSE]
        for (ii in 1:nchecks) {
          if (any(duplicated(centers[[length(centers)]]))) {
            centers[[length(centers)]] <- x[sample(1:n, k),,drop=FALSE]
          } else break;
        }
        if (any(duplicated(centers[[length(centers)]]))) {
          stop("Having trouble finding non-duplicated centers.\n")
        }
      }
    }
    return(centers)
  }

  ####################################################################
  # A function for aggregating results as foreach is running, to avoid
  # memory overhead.
  choosebest <- function(a, b) {
    if (!is.na(sum(a$withinss)) & is.na(sum(b$withinss))) return(a)
    else if (is.na(sum(a$withinss)) & !is.na(sum(b$withinss))) return(b)
    if ( sum(a$withinss) < sum(b$withinss) ) {
      return(a)
    } else {
      return(b)
    }
  }
    
  #################################################
  # Check centers for sanity and consider nstart>1:
  if (!is.matrix(centers)) {
    if (is.numeric(centers) && length(centers)==1 && centers>0) {
      k <- centers
      centers <- getcenters(x, k, nstart)
    } else stop("centers must be a matrix of centers or number of clusters > 0")
  } else {
    k <- nrow(centers)
    if (nstart>1) {
      warning(paste("Random starting points will be used",
                    "(not the centers you provided), because nstart>1.\n"))
      centers <- getcenters(x, k, nstart)
    } else {
      if (any(duplicated(centers))) {
        stop(paste("Error: if you provide centers,",
                   "they had better not have duplicates.\n"))
      }
      centers <- list(centers)
    }
  }

  ###############################################################
  # At this point, centers is a list of length nstart of matrices
  # of starting centers without duplicates.
  # I think I allow k=1 cluster, too, but check it later.
  # Note that if number of columns is HUGE, the centers will
  # be memory-intensive.

  if (is.matrix(x)) {
    if (getDoParWorkers()>1) {
      # Generate big.matrix copy so we can work in parallel; we
      # assume in-memory is fine given that x is a matrix.
      if (is.integer(x)) {
        y <- as.big.matrix(x, type="integer")
      } else {
        y <- as.big.matrix(x, type="double")
      }
      xdesc <- describe(y)
    } else { xdesc <- NULL }       # This is the signal of a matrix.
  } else {
    if (is.big.matrix(x)) {
      xdesc <- describe(x)
    } else {
      stop("x must be a matrix or a big.matrix.\n")
    }
  }

  nr <- nrow(x)
  if (typeof(x)=="char") mattype <- 1
  if (typeof(x)=="short") mattype <- 2
  if (typeof(x)=="integer") mattype <- 4
  if (typeof(x)=="double") mattype <- 8

  cen <- NA # Note this is because of foreach's interaction with the
            # R check parser.  There is no real problem, but I need this
            # to avoid a warning.

  # Do the work, possibly in parallel with nstart>1 and a registered
  # parallel backend.
  ans <- foreach(cen=centers, .combine="choosebest", 
                 .packages=c("bigmemory", "biganalytics")) %dopar% {
    center <- big.matrix(nrow(cen), ncol(cen), type="double")
    center[,] <- cen
    clust <- big.matrix(nr, 1, type="integer")
    clustsizes <- big.matrix(nrow(cen), 1, type="double")
    wss <- big.matrix(nrow(cen), 1, type="double")
    if (is.null(xdesc)) {
      # .Call with the matrix, which has to be either integer or double.
      if (mattype==4) {
        res <- .Call("kmeansRIntMatrix", x,
                     center@address, clust@address, clustsizes@address,
                     wss@address, as.integer(iter.max), as.integer(dist_calc),
                     PACKAGE="biganalytics")
      } else {
        res <- .Call("kmeansRNumericMatrix", x,
                     center@address, clust@address, clustsizes@address,
                     wss@address, as.integer(iter.max), as.integer(dist_calc),
                     PACKAGE="biganalytics")
      }
    } else {
      # .Call with the big.matrix
      x <- attach.big.matrix(xdesc)
      res <- .Call("kmeansBigMatrix", x@address,
                   center@address, clust@address, clustsizes@address,
                   wss@address, as.integer(iter.max), as.integer(dist_calc),
                   PACKAGE="biganalytics")
    }

    temp <- list(cluster=clust[,],
                 centers=center[,],
                 withinss=wss[,],
                 size=clustsizes[,],
                 iters=res)

    return(temp)

  } # End of the foreach() body.

  if (ans$iters>=iter.max) 
    warning("bigkmeans did not converge in ", iter.max, " iterations.\n")
  ans$iters <- NULL                     # This removes this from the list.
  class(ans) <- "kmeans"

  return(ans)

}


