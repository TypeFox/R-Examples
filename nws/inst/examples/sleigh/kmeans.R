library(nws)

kmeansInit <-
function()
{
  # fetch our chunk of the data, and put it in the global variable "x"
  x <<- nwsFetch(SleighUserNws, paste("x", SleighRank + 1, sep=""))
  invisible()
}

kmeansWorker <-
function(cen, algorithm)
{
  # x is a global variable that was initialized by kmeansInit
  km <- kmeans(x, cen, iter.max=1, algorithm=algorithm)
  km$id <- SleighRank + 1  # sneak in the id of the worker
  km
}

kmeansNws <-
function(s, x, centers, iter.max=10,
         algorithm=c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
{
  n <- workerCount(s)
  x <- as.matrix(x)
  if(missing(centers))
    stop("'centers' must be a matrix")
  m <- nrow(x)
  k <- nrow(centers)
  algorithm <- match.arg(algorithm)

  # send matrix to nws server in "n" chunks
  for (proc in seq(length.out=n))
    nwsStore(s@userNws, paste("x", proc, sep=""),
             x[seq(proc, m, n), , drop=FALSE])

  # each worker fetches its own chunk of "x" from the Sleigh workspace
  eachWorker(s, kmeansInit)

  cluster <- integer(m)
  for (i in seq(length.out=iter.max)) {
    # get the results for the next iteration from each chunk of data
    results <- eachWorker(s, kmeansWorker, centers, algorithm)

    # compute centers, size, and cluster from the partial results
    tmpcenters <- matrix(0, k, ncol(centers))
    tmpsize <- double(k)
    tmpcluster <- integer(m)

    for (cl in results) {
      tmp <- cl$centers * cl$size
      tmp[cl$size == 0,] <- 0  # zero out rows corresponding to empty clusters
      tmpcenters <- tmpcenters + tmp
      tmpsize <- tmpsize + cl$size
      tmpcluster[seq(cl$id, m, n)] <- cl$cluster
    }

    # check for convergence
    if (all(cluster == tmpcluster))
      break

    centers <- tmpcenters / tmpsize  # returns NaNs for empty clusters
    cluster <- tmpcluster

    # if this is the last iteration, then we failed to converge
    if (i == iter.max) {
      warning("did not converge in ", i, " iterations", call.=FALSE)
      i <- i + 1  # this indicates that it didn't converge
    }
  }

  # withinss computation is done here to avoid doing it more than once
  withinss <- rep(0, k)
  for (cl in results)
    withinss <- withinss + cl$withinss

  # generate the "kmeans" object (with an extra "niter" tag)
  out <- list(cluster=cluster, centers=centers, withinss=withinss,
              size=tmpsize, niter=i)
  class(out) <- "kmeans"
  out
}

kmeansTest <-
function(s)
{
  if (missing(s))
    error("kmeansTest needs a sleigh object")

  numcol <- 2
  k <- 4
  x <- rbind(matrix(rnorm(40, mean=0, sd=0.6), ncol=numcol),
             matrix(rnorm(40, mean=1, sd=0.6), ncol=numcol))

  iter.max <- 10
  colnames(x) <- 1:numcol

  # pick starting point for cluster centers randomly
  cn <- unique(x)
  centers <- cn[sample(1:nrow(cn), k), , drop=FALSE]
  rownames(centers) <- 1:k

  cl <- kmeansNws(s, x, centers, iter.max=iter.max, algorithm="Lloyd")

  # display the results
  print(cl)
  plot(x, col=cl$cluster)
  points(cl$centers, col=1:k, pch=8, cex=2)

  # compare to sequential results
  scl <- kmeans(x, centers, iter.max=iter.max, algorithm="Lloyd")
  print(cl$centers)
  print(scl$centers)
  if (!identical(all.equal(cl$centers, scl$centers), TRUE))
    stop("parallel centers do not match sequential centers")
  cat("WSS comparison: ")
  print(all.equal(cl$withinss, scl$withinss))
}

s <- sleigh()
kmeansTest(s)
