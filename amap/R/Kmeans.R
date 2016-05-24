Kmeans <-
function(x, centers, iter.max = 10, nstart = 1,
         method = "euclidean")
{
  dokmeans <- function()
    {
      Z <- .C("kmeans_Lloyd2", as.double(x), as.integer(m),
              as.integer(ncol(x)),
              centers = as.double(centers), as.integer(k),
              c1 = integer(m), iter = as.integer(iter.max),
              nc = integer(k), wss = double(k),
              method=as.integer(method),
              PACKAGE="amap")
      if (Z$iter > iter.max) 
        warning("did not converge in ", iter.max, " iterations", 
                call. = FALSE)
      if (any(Z$nc == 0)) 
        warning("empty cluster: try a better set of initial centers", 
                call. = FALSE)
      Z
    }


  
  METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
               "binary","pearson","correlation","spearman","kendall","abspearson","abscorrelation")
  method <- pmatch(method, METHODS)
  if (is.na(method)) 
    stop("invalid distance method")
  if (method == -1) 
    stop("ambiguous distance method")
  
  if(class(x) == "exprSet")
  { 
    library(Biobase)
     x <- Biobase::exprs(x)
  }

  x <- as.matrix(x)
  m <- nrow(x)
  if(missing(centers))
    stop("'centers' must be a number or a matrix")
  if(length(centers) == 1) {
    k <- centers
    ## we need to avoid duplicates here
    if(nstart == 1)
      centers <- x[sample(1 : m, k), , drop = FALSE]
    if(nstart >= 2 || any(duplicated(centers))) {
      cn <- unique(x)
      mm <- nrow(cn)
      if(mm < k)
        stop("more cluster centers than distinct data points.")
      centers <- cn[sample(1:mm, k), , drop=FALSE]
    }
  } else {
    centers <- as.matrix(centers)
    if(any(duplicated(centers)))
      stop("initial centers are not distinct")
    cn <- NULL
    k <- nrow(centers)
    if(m < k)
      stop("more cluster centers than data points")
  }
  if(iter.max < 1) stop("'iter.max' must be positive")
  if(ncol(x) != ncol(centers))
    stop("must have same number of columns in 'x' and 'centers'")
  
  
  Z <- .C("kmeans_Lloyd2", as.double(x), as.integer(m),
          as.integer(ncol(x)),
          centers = as.double(centers), as.integer(k),
          c1 = integer(m), iter = as.integer(iter.max),
          nc = integer(k), wss = double(k),
          method=as.integer(method),
          PACKAGE="amap")
  if(Z$iter > iter.max)
    warning("did not converge in ",
            iter.max, " iterations", call.=FALSE)
  if(any(Z$nc == 0))
    warning("empty cluster: try a better set of initial centers", call.=FALSE)
    
  if(nstart >= 2 && !is.null(cn)) {
    best <- sum(Z$wss)
    for(i in 2:nstart) {
      centers <- cn[sample(1:mm, k), , drop=FALSE]
      ZZ <- dokmeans()
      if((z <- sum(ZZ$wss)) < best) {
        Z <- ZZ
        best <- z
      }
    }
  }
  centers <- matrix(Z$centers, k)
  dimnames(centers) <- list(1:k, dimnames(x)[[2]])
  cluster <- Z$c1
  if(!is.null(rn <- rownames(x)))
    names(cluster) <- rn
  
  out <- list(cluster = cluster, centers = centers, withinss = Z$wss,
                size = Z$nc)
  class(out) <- "kmeans"
  out
}

