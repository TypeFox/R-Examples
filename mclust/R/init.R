#############################################################################
## Initialization for d-dim data ############################################
#############################################################################

# This new version allowing transformation of the data. 
# By default it behaves as the old function
hc <- function(data, modelName = mclust.options("hcModelNames")[1], 
               use = mclust.options("hcUse"), ...)
{
  if(!any(modelName == c("E", "V", "EII", "VII", "EEE", "VVV")))
     stop("invalid 'modelName' argument for model-based hierarchical clustering")

  if(!any(use == c("VARS", "STD", "SPH", "PCS", "PCR", "SVD")))
     stop("invalid 'use' argument for  model-based hierarchical clustering")

  funcName <- paste("hc", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc$use <- mc$modelName <- NULL
  data <- data.matrix(data)

  dropCols <- function(x)
  { # select only those columns of matrix x with all finite numeric values
    x[,apply(x, 2, function(x) all(is.finite(x))), drop = FALSE]
  }

  use <- toupper(use)
  switch(use,
         "VARS" = { Z <- data },
         "STD" = { Z <- scale(data, center = TRUE, scale = TRUE) 
                   Z <- dropCols(Z) },
         "PCR" = { data <- scale(data, center = TRUE, scale = TRUE)
                   data <- dropCols(data)
                   SVD <- svd(data, nu=0)
                   # evalues <- sqrt(SVD$d^2/(nrow(data)-1))
                   Z <- data %*% SVD$v },
         "PCS" = { data <- scale(data, center = TRUE, scale = FALSE)
                   SVD <- svd(data, nu=0)
                   # evalues <- sqrt(SVD$d^2/(nrow(data)-1))
                   Z <- data %*% SVD$v 
                   Z <- dropCols(Z) },
         "SPH" = { data <- scale(data, center = TRUE, scale = FALSE)
                   n <- nrow(data); p <- ncol(data)
                   Sigma <- var(data) * (n - 1)/n
                   SVD <- svd(Sigma, nu = 0)
                   Z <- data %*% SVD$v %*% diag(1/sqrt(SVD$d), p, p) 
                   Z <- dropCols(Z) },
         "SVD" = { data <- scale(data, center = TRUE, scale = TRUE)
                   data <- dropCols(data)
                   n <- nrow(data); p <- ncol(data)
                   SVD <- svd(data, nu=0)
                   Z <- data %*% SVD$v %*% diag(1/sqrt(SVD$d), p, p) },
         stop("'use' argument not allowed. See help(mclust.options)")
            )
  # call the proper hc<funcName> function
  mc$data <- Z 
  mc[[1]] <- as.name(funcName)
  out <- eval(mc, parent.frame())
  attr(out, "call") <- match.call()
  class(out) <- "hc"
  return(out)
}


print.hc <- function(x, ...) 
{
  if(!is.null(attr(x, "call"))) 
    cat("Call:\n", deparse(attr(x, "call")), "\n\n", sep = "")
  cat("Model-Based Agglomerative Hierarchical Clustering:\n")
  if(!is.null(attr(x, "modelName")))
    cat("Model name        = ", attr(x, "modelName"), "\n")
  if(!is.null(attr(x, "dimensions")))
    cat("Number of objects = ", attr(x, "dimensions")[1], "\n")
  invisible(x)
}

randomPairs <- function(data, seed, ...)
{
  if(!missing(seed)) set.seed(seed)
  data <- as.matrix(data)
  n <- nrow(data)
  tree <- matrix(sample(1:n, n, replace = FALSE), nrow = 2, ncol = ceiling(n/2))
  tree <- apply(tree, 2, sort)
  ind <- unique(tree[1,])
  while(ncol(tree) < (n-1))
  { 
    addtree <- sort(sample(ind, size = 2, replace = FALSE))
    ind <- setdiff(ind, addtree[2])
    tree <- cbind(tree, addtree)
  }
  dimnames(tree) <- NULL
  structure(tree, initialPartition = 1:n, dimensions = c(n,2))
}

hclass <- function(hcPairs, G)
{
  initial <- attributes(hcPairs)$init
  n <- length(initial)
  k <- length(unique(initial))
  G <- if(missing(G)) k:2 else rev(sort(unique(G)))
  select <- k - G
  if(length(select) == 1 && !select)
    return(matrix(initial, ncol = 1, dimnames = list(NULL, 
                                                     as.character(G))))
  bad <- select < 0 | select >= k
  if(all(bad))
    stop("No classification with the specified number of clusters")
  if(any(bad) & mclust.options("warn"))
    { warning("Some selected classifications are inconsistent with mclust object") }
  L <- length(select)
  cl <- matrix(as.double(NA), nrow = n, ncol = L, 
               dimnames = list(NULL, as.character(G)))
  if(select[1])
    m <- 1
  else {
    cl[, 1] <- initial
    m <- 2
  }
  for(l in 1:max(select)) {
    ij <- hcPairs[, l]
    i <- min(ij)
    j <- max(ij)
    initial[initial == j] <- i
    if(select[m] == l) {
      cl[, m] <- initial
      m <- m + 1
    }
  }
  apply(cl[, L:1, drop = FALSE], 2, partconv, consec = TRUE)
}

hcEII <- function(data, partition, minclus = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #====================================================================
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || length(dimdat[dimdat > 1]) == 1)
  #if(oneD || length(dimdat) > 2)
  #  stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    { stop("initial number of clusters is not greater than minclus") }
  if(n <= p & mclust.options("warn"))
    { warning("# of observations <= data dimension") }
  #=============================================================
  storage.mode(data) <- "double"
  ld <- max(c((l * (l - 1))/2, 3 * m))
  temp <- .Fortran("hceii",
                   data,
                   as.integer(n),
                   as.integer(p),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   double(p),
                   as.integer(ld),
                   double(ld),
                   PACKAGE = "mclust")[c(1, 9)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
  temp[[2]] <- temp[[2]][1:m]
  structure(t(temp[[1]]), initialPartition = partition, 
            dimensions = dimdat, modelName = "EII", 
            call =  match.call())
}

hcEEE <- function(data, partition, minclus = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #=====================================================================
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || length(dimdat[dimdat > 1]) == 1)
  #if(oneD || length(dimdat) > 2)
  #  stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(n <= p & mclust.options("warn"))
    warning("# of observations <= data dimension")
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  storage.mode(data) <- "double"
  
  ## R 2.12.0: 32 bit Windows build fails due to compiler bug
  ## workaround: removal (hopefully temporary) of hc functionality for EEE
  # Luca: commented the next line and uncommented below
  #  stop("hc for EEE model is not currently supported")
  
  temp <- .Fortran("hceee",
                   data,
                   as.integer(n),
                   as.integer(p),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   if(p < 3) integer(m) else integer(1),
                   if(p < 4) integer(m) else integer(1),
                   double(p),
                   double(p * p),
                   double(p * p),
                   double(p * p),
                   PACKAGE = "mclust")[c(1, 7:10)]
  #
  # currently temp[[5]] is not output
  temp[[4]] <- temp[[4]][1:2]
  temp[[5]] <- temp[[5]][1:2]
  names(temp[[5]]) <- c("determinant", "trace")
  temp[[1]] <- temp[[1]][1:(m + 1),  ]
  if(p < 3)
    tree <- rbind(temp[[2]], temp[[3]])
  else if(p < 4)
    tree <- rbind(temp[[1]][-1, 3], temp[[3]])
  else tree <- t(temp[[1]][-1, 3:4, drop = FALSE])
  determinant <- temp[[1]][, 1]
  attr(determinant, "breakpoints") <- temp[[4]]
  trace <- temp[[1]][, 2]
  structure(tree,  initialPartition = partition, 
            dimensions = dimdat, modelName = "EEE", 
            call = match.call())
}

hcVII <- function(data, partition, minclus = 1, alpha = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #=====================================================================
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || length(dimdat[dimdat > 1]) == 1)
  #if(oneD || length(dimdat) > 2)
  #  stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(n <= p & mclust.options("warn"))
    warning("# of observations <= data dimension")
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  storage.mode(data) <- "double"
  ll <- (l * (l - 1))/2
  ld <- max(n, ll, 3 * m)
  alpha <- alpha * traceW(data/sqrt(n * p))
  alpha <- max(alpha, .Machine$double.eps)
  temp <- .Fortran("hcvii",
                   data,
                   as.integer(n),
                   as.integer(p),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   as.double(alpha),
                   double(p),
                   as.integer(ld),
                   double(ld),
                   PACKAGE = "mclust")[c(1, 10)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
  temp[[2]] <- temp[[2]][1:m]
  structure(t(temp[[1]]), initialPartition = partition, 
            dimensions = dimdat, modelName = "VII", 
            call = match.call())
}

hcVVV <- function(data, partition, minclus = 1, alpha = 1, beta = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || length(dimdat[dimdat > 1]) == 1)
  #if(oneD || length(dimdat) > 2)
  #  stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(n <= p & mclust.options("warn"))
    warning("# of observations <= data dimension")
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  storage.mode(data) <- "double"
  ll <- (l * (l - 1))/2
  #  dp <- duplicated(partition)
  #x[c((1:n)[!dp],(1:n)[dp]),], 
  #as.integer(c(partition[!dp], partition[dp])), 
  ld <- max(n, ll + 1, 3 * m)
  alpha <- alpha * traceW(data/sqrt(n * p))
  alpha <- max(alpha, .Machine$double.eps)
  temp <- .Fortran("hcvvv",
                   cbind(data, 0.),
                   as.integer(n),
                   as.integer(p),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   as.double(alpha),
                   as.double(beta),
                   double(p),
                   double(p * p),
                   double(p * p),
                   double(p * p),
                   as.integer(ld),
                   double(ld),
                   PACKAGE = "mclust")[c(1, 14)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
  temp[[2]] <- temp[[2]][1:m]
  structure(t(temp[[1]]), initialPartition = partition, 
            dimensions = dimdat, modelName = "VVV", 
            call = match.call())
}

## Initialization for 1-dim data ############################################

# This version is bugged when a quantile is equal to the following
# qclass <- function (x, k)
# {
#   q <- quantile(x, seq(from = 0, to = 1, by = 1/k))
#   cl <- rep(0, length(x))
#   q[1] <- q[1] - 1
#   for(i in 1:k)
#     cl[x > q[i] & x <= q[i+1]] <- i
#   return(cl)
# }

# This should correct the above bug
qclass <- function (x, k) 
{
  x <- as.vector(x)
  # eps <- sqrt(.Machine$double.eps) 
  # numerical accuracy problem if scale of x is large, so make tolerance
  # scale dependent
  eps <- sd(x)*sqrt(.Machine$double.eps)
  q <- NA
  n <- k
  while(length(q) < (k+1))
  { n <- n + 1
    q <- unique(quantile(x, seq(from = 0, to = 1, length = n))) 
  }
  if(length(q) > (k+1))
  { dq <- diff(q)
    nr <- length(q)-k-1
    q <- q[-order(dq)[1:nr]]
  }
  q[1] <- min(x) - eps
  q[length(q)] <- max(x) + eps
  cl <- rep(0, length(x))
  for(i in 1:k) 
     { cl[ x >= q[i] & x < q[i+1] ] <- i }
  return(cl)
}

hcE <- function(data, partition, minclus = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #====================================================================
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  storage.mode(data) <- "double"
  ld <- max(c((l * (l - 1))/2, 3 * m))
  temp <- .Fortran("hc1e",
                   data,
                   as.integer(n),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   as.integer(ld),
                   double(ld),
                   PACKAGE = "mclust")[c(1, 3, 7)]
  temp[[1]] <- temp[[1]][1:m]
  temp[[2]] <- temp[[2]][1:m]
  temp[[3]] <- temp[[3]][1:m]
  structure(rbind(temp[[1]], temp[[2]]),   initialPartition = partition, 
            dimensions = n, modelName = "E",
            call = match.call())
}

hcV <- function(data, partition, minclus = 1, alpha = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #=====================================================================
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  storage.mode(data) <- "double"
  alpha <- alpha * (vecnorm(data - mean(data))^2/n)
  alpha <- min(alpha, .Machine$double.eps)
  ld <- max(c((l * (l - 1))/2, 3 * m))
  temp <- .Fortran("hc1v",
                   data,
                   as.integer(n),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   as.double(alpha),
                   as.integer(ld),
                   double(ld),
                   PACKAGE = "mclust")[c(1, 3, 8)]
  temp[[1]] <- temp[[1]][1:m]
  temp[[2]] <- temp[[2]][1:m]
  temp[[3]] <- temp[[3]][1:m]
  structure(rbind(temp[[1]], temp[[2]]),   initialPartition = partition, 
            dimensions = n, modelName = "V",
            call = match.call())
}
