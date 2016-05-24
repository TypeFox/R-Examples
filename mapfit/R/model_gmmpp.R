
gmmpp <- function(size, alpha, D0, D1, class="dgeMatrix") {
  if (missing(size)) {
    if (missing(alpha) || missing(D0) || missing(D1)) {
      stop("alpha, D0 and D1 are needed.")
    } else {
      size <- length(alpha)
    }
  } else {
    if (!missing(alpha) && !missing(D0) && !missing(D1)) {
      warning("size is ignored.")
      size <- length(alpha)
    } else {
      if (missing(alpha)) {
        # warning("alpha is set by default")
        alpha <- rep(1.0/size, size)
      }
      if (missing(D1)) {
        # warning("D1 is set by default")
        D1 <- diag(1.0, size, size)
      }
      if (missing(D0)) {
        # warning("D0 is set by default")
        D0 <- matrix(1.0, size, size)
        diag(D0) <- -(rep(size-1, size) + apply(D1,1,sum))
      }
    }
  }
  if (!is(D0, class)) {
    D0 <- as(D0, class)
  }
  if (!is(D1, class)) {
    D1 <- as(D1, class)
  }
  zero <- 1.0e-8
  df <- sum(abs(alpha) > zero) - 1 + sum(abs(D0) > zero) +
    sum(abs(D1) > zero) - size + sum(abs(diag(D0)) < zero)
  new("gmmpp", size=size, alpha=alpha, D0=D0, D1=D1, df=df)
}

gmmpp.tridiag <- function(size, class="dgeMatrix") {
  if (size <= 2) {
    gmmpp(size, class=class)
  } else {
    alpha <- rep(1.0/size, size)
    D0 <- matrix(0, size, size)
    D1 <- diag(1, size)
    D0[1,1] <- -2
    D0[1,2] <- 1
    for (i in 2:(size-1)) {
      D0[i,i] <- -3
      D0[i,i-1] <- 1
      D0[i,i+1] <- 1
    }
    D0[size,size] <- -2
    D0[size,size-1] <- 1
    gmmpp(alpha=alpha, D0=D0, D1=D1, class=class)
  }
}

## S4 methods

setMethod("emfit.init", signature(model = "gmmpp", data = "mapdata"),
  function(model, data, verbose = list(), ...) {
    gmmpp.param.kmeans(size=model@size, data=data@data, 
      skelal=model@alpha, skelD0=as.matrix(model@D0),
      skelD1=as.matrix(model@D1), verbose=verbose, ...)
  }
)

gmmpp.param.kmeans <- function(size, data, skelal, skelD0, skelD1, verbose, ...) {
  if (missing(skelal)) skelal <- rep(1, size)
  if (missing(skelD0)) skelD0 <- matrix(1, size, size)
  if (missing(skelD1)) skelD1 <- matrix(1, size, size)

  maxt <- max(data$time)
  maxg <- max(data$counts + data$instant)
  x <- cbind(data$time, (data$counts + data$instant))
  v <- cbind(data$time / maxt, (data$counts + data$instant) / maxg)
  result <- kmeans(v, size)
  
  diagelem <- sapply(1:size, function(k) {
    (sum(x[result$cluster == k,2]) + 1) /
    sum(x[result$cluster == k,1])})

  alpha <- skelal * runif(size)
  alpha <- alpha / sum(alpha)

  diag(skelD0) <- 0
  D0 <- skelD0 * matrix(runif(size*size), size, size)
  D1 <- skelD1 * matrix(runif(size*size), size, size)
  d <- diagelem / (apply(D0, 1, sum) + apply(D1, 1, sum))
  D0 <- D0 * d
  D1 <- D1 * d
  diag(D0) <- -diagelem
  gmmpp(alpha=alpha, D0=D0, D1=D1, ...)
}

#### estep

setMethod("emfit.estep", signature(model = "gmmpp", data = "mapdata"),
  function(model, data, eps = 1.0e-8, devide = 20, ...) {
    res <- .Call(mapfit_estep_gmmpp, model, data, eps, devide)
    list(eres=list(eb=res[[1]], ez=res[[2]], en0=res[[3]], en1=res[[4]]), llf=res[[5]])
  })

