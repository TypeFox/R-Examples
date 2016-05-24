cclust <- function (x, centers, iter.max = 100, verbose = FALSE, dist = "euclidean", 
    method = "kmeans", rate.method = "polynomial", rate.par = NULL) 
{
    xrows <- dim(x)[1]
    xcols <- dim(x)[2]
    xold <- x
    perm <- sample(xrows)
    x <- x[perm, ]
    # initial values are given
    if (is.matrix(centers)) 
        ncenters <- dim(centers)[1]
    else {
        # take centers random vectors as initial values
        ncenters <- centers
        centers <- x[rank(runif(xrows))[1:ncenters], ]
    }
    dist <- pmatch(dist, c("euclidean", "manhattan"))
    if (is.na(dist)) 
        stop("invalid distance")
    if (dist == -1) 
      stop("ambiguous distance")
    method <- pmatch(method, c("kmeans", "hardcl", "neuralgas"))
    if (is.na(method)) 
        stop("invalid clustering method")
    if (method == -1) 
        stop("ambiguous clustering method")
    rate.method <- pmatch(rate.method, c("polynomial", "exponentially.decaying"))
    if (is.na(rate.method)) 
        stop("invalid learning rate method")
    if (rate.method == -1) 
        stop("ambiguous learning rate method")
    if (method == 2) {
        if (rate.method == 1 && missing(rate.par)) {
            rate.par <- c(1e-00, 0e-00)
        }
        else if (rate.method == 2 && missing(rate.par)) {
            rate.par <- c(0.1, 1e-04)
        }
    }
    if (method == 3 && missing(rate.par)) {
        rate.par <- c(0.5, 0.005, 10, 0.01)
    }
    initcenters <- centers
   # dist <- matrix(0, xrows, ncenters)
    # necessary for empty clusters
    pos <- as.factor(1:ncenters)
    rownames(centers) <- pos
    iter <- integer(1)
    changes <- integer(iter.max)
    cluster <- integer(xrows)
    clustersize <- integer(ncenters)
    if (method == 1) {
      retval <- .C("kmeans", xrows = as.integer(xrows),
                   xcols = as.integer(xcols), 
                   x = as.double(x), ncenters = as.integer(ncenters), 
                   centers = as.double(centers),
                   cluster = as.integer(cluster), 
                   iter.max = as.integer(iter.max), iter = as.integer(iter), 
                   changes = as.integer(changes),
                   clustersize = as.integer(clustersize), 
                   verbose = as.integer(verbose),
                   dist = as.integer(dist-1), PACKAGE="cclust")
    }
    else if (method == 2) {
      retval <- .C("hardcl", xrows = as.integer(xrows), xcols = as.integer(xcols), 
                   x = as.double(x), ncenters = as.integer(ncenters), 
                   centers = as.double(centers),
                   cluster = as.integer(cluster), 
                   iter.max = as.integer(iter.max), iter = as.integer(iter), 
                   clustersize = as.integer(clustersize),
                   verbose = as.integer(verbose), 
                   dist = as.integer(dist-1),
                   methrate = as.integer(rate.method-1), 
                   par = as.double(rate.par), PACKAGE="cclust")
    }
    else if (method == 3) {
      retval <- .C("neuralgas", xrows = as.integer(xrows), 
                   xcols = as.integer(xcols), x = as.double(x),
                   ncenters = as.integer(ncenters), 
                   centers = as.double(centers),
                   cluster = as.integer(cluster), 
                   iter.max = as.integer(iter.max), iter = as.integer(iter), 
                   clustersize = as.integer(clustersize),
                   verbose = as.integer(verbose), 
                   dist = as.integer(dist-1), par = as.double(rate.par),
                   PACKAGE="cclust")
    }
    centers <- matrix(retval$centers, ncol = xcols, dimnames = dimnames(initcenters))
    cluster <- retval$cluster + 1
    cluster <- cluster[order(perm)]
    if (method == 1) {
        methrate <- NA
        par <- NA
    }
  if (method == 3) {
        methrate <- NA
    }

    withinss <- function(clobj, x){
      
      retval <- rep(0, nrow(clobj$centers))
      x <- (x - clobj$centers[clobj$cluster, ])^2
      for(k in 1:nrow(clobj$centers)){
        retval[k] <- sum(x[clobj$cluster==k,])
      }
      retval
    }
    
    within <- withinss(list(centers = centers, cluster = cluster), xold)
    
    retval <- list(centers = centers, initcenters = initcenters, 
        ncenters = ncenters, cluster = cluster, size = retval$clustersize, 
        iter = retval$iter - 1, changes = retval$changes, dist = dist, 
        method = method, rate.method = rate.method, rate.par = rate.par, 
        call = match.call(), withinss = within)
    class(retval) <- c("cclust")
    return(retval)
  }



print.cclust <- function (x, ...)
  {
    clobj <- x
    if (!is.null(clobj$iter))
      cat("\n                            Clustering on Training Set\n\n\n")
    else
      cat("\n                              Clustering on Test Set\n\n\n")
    
    cat("Number of Clusters: ", clobj$ncenters, "\n")
    cat("Sizes  of Clusters: ", clobj$size, "\n\n")
    if (clobj$method!=1)
      cat("Learning Parameters:",clobj$rate.par,"\n\n")
    
  if (clobj$method==1)
    {
  if (!is.null(clobj$iter))
      {
        if (clobj$iter < length(clobj$changes))
          cat("Algorithm converged after", clobj$iter, "iterations.\n")
        else
          cat("Algorithm did not converge after", clobj$iter, "iterations.\n")
        cat("Changes:", clobj$changes[1:clobj$iter], "\n\n")
      }
    }
 
  }


predict.cclust <- function(object, newdata, ...){

  clobj <- object
  x <- newdata
  xrows<-dim(x)[1]
  xcols<-dim(x)[2]
  ncenters <- clobj$ncenters
  cluster <- integer(xrows)
  clustersize <- integer(ncenters)
  

  if(dim(clobj$centers)[2] != xcols){
    stop("Number of variables in cluster object and x are not the same!")
  }

  
  retval <- .C("assign",
               xrows = as.integer(xrows),
               xcols = as.integer(xcols),
               x = as.double(x),
               ncenters = as.integer(ncenters),
               centers = as.double(clobj$centers),
               cluster = as.integer(cluster),
               clustersize = as.integer(clustersize),
               dist = as.integer(clobj$dist-1),PACKAGE="cclust")

  
     

  clobj$initcenters <- NULL
  clobj$iter <- NULL
  clobj$changes <- NULL
  clobj$cluster <- retval$cluster+1
  clobj$size <- retval$clustersize

  return(clobj)
}



