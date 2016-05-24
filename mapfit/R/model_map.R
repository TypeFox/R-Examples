
map <- function(size, alpha, D0, D1, class="CsparseMatrix") {
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
        D1 <- matrix(1.0, size, size)
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
  new("map", size=size, alpha=alpha, D0=D0, D1=D1, df=df)
}

mmpp <- function(size, class="CsparseMatrix") {
  alpha <- rep(1.0/size, size)
  D1 <- diag(1.0, size)
  D0 <- matrix(1.0, size, size)
  diag(D0) <- -(rep(size-1, size) + apply(D1,1,sum))
  map(alpha=alpha, D0=D0, D1=D1, class=class)
}

map.bidiag <- function(size, class="CsparseMatrix") {
  if (size <= 1) {
    map(size, class=class)
  } else {
    alpha <- rep(1.0/size, size)
    D0 <- matrix(0, size, size)
    D1 <- matrix(1, size, size)
    for (i in 1:(size-1)) {
      D0[i,i] <- -1-size
      D0[i,i+1] <- 1
    }
    D0[size,size] <- -size
    map(alpha=alpha, D0=D0, D1=D1, class=class)
  }
}

map.tridiag <- function(size, class="CsparseMatrix") {
  if (size <= 2) {
    map(size, class=class)
  } else {
    alpha <- rep(1.0/size, size)
    D0 <- matrix(0, size, size)
    D1 <- matrix(1, size, size)
    D0[1,1] <- -1-size
    D0[1,2] <- 1
    for (i in 2:(size-1)) {
      D0[i,i] <- -2-size
      D0[i,i-1] <- 1
      D0[i,i+1] <- 1
    }
    D0[size,size] <- -1-size
    D0[size,size-1] <- 1
    map(alpha=alpha, D0=D0, D1=D1, class=class)
  }
}

mmpp.tridiag <- function(size, class="CsparseMatrix") {
  if (size <= 2) {
    map(size, class=class)
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
    map(alpha=alpha, D0=D0, D1=D1, class=class)
  }
}

map.mmoment <- function(k, map) {
  D0 <- as.matrix(map@D0)
  D1 <- as.matrix(map@D1)
  piv <- as.vector(ctmc.st(D0+D1) %*% D1)
  piv <- piv / sum(piv)
  tmp <- piv
  tmp2 <- 1.0
  res <- numeric(0)
  for (i in 1:k) {
    tmp <- msolve(alpha=1.0, A=-D0, x=tmp, transpose=TRUE)
    tmp2 <- tmp2 * i
    res <- c(res, tmp2 * sum(tmp))
  }
  res
}

map.jmoment <- function(lag, map) {
  D0 <- as.matrix(map@D0)
  D1 <- as.matrix(map@D1)
  piv <- as.vector(ctmc.st(D0+D1) %*% D1)
  piv <- piv / sum(piv)
  P <- mpow(mpow(A=-D0, m=-1) %*% D1, m=lag)
  vone <- rep(1, map@size)

  fmat <- matrix(0, map@size, map@size)
  fmat[,1] <- piv
  for (i in 2:map@size) {
    fmat[,i] <- msolve(alpha=1.0, A=-D0, x=fmat[,i-1], transpose=TRUE)
  }

  res <- matrix(0, map@size, map@size)
  tmp <- vone
  for (j in 1:map@size) {
    for (i in 1:map@size) {
      res[i,j] <- gamma(i) * gamma(j) * as.vector(fmat[,i] %*% P %*% tmp)
    }
    tmp <- msolve(alpha=1.0, A=-D0, x=tmp)
  }
  res
}

map.acf <- function(map) {
  D0 <- as.matrix(map@D0)
  D1 <- as.matrix(map@D1)
  piv <- as.vector(ctmc.st(D0+D1) %*% D1)
  piv <- piv / sum(piv)
  P <- mpow(A=-D0, m=-1) %*% D1
  vone <- rep(1, map@size)

  piv <- msolve(alpha=1.0, A=-D0, x=piv, transpose=TRUE)
  vone <- msolve(alpha=1.0, A=-D0, x=vone)

  mres <- map.mmoment(2, map)
  (sapply(1:map@size, function(k) as.vector(piv %*% mpow(P,k) %*% vone)) - mres[1]^2) / (mres[2] - mres[1]^2)
}

## S4 methods

setMethod("emfit.print", signature(model = "map"),
  function(model, ...) {
    cat(gettextf("Size : %d\n", model@size))
    cat("Initial : ", model@alpha, "\n")
    cat("Infinitesimal generator : \n")
    print(model@D0)
    print(model@D1)
  }
)

setMethod("emfit.df", signature(model = "map"),
  function(model, stationary = FALSE, ...) {
    if (stationary)
      model@df - model@size + 1
    else
      model@df
  }
)

## init

setMethod("emfit.init", signature(model = "map", data = "mapdata"),
  function(model, data, verbose = list(), ...) {
    map.param.kmeans(size=model@size, data=data@data, 
      skelal=model@alpha, skelD0=as.matrix(model@D0),
      skelD1=as.matrix(model@D1), verbose=verbose, class=class(model@D0), ...)
  }
)

map.param.kmeans <- function(size, data, skelal, skelD0, skelD1, verbose, class, ...) {
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
  map(alpha=alpha, D0=D0, D1=D1, class=class)
}

#### estep

setMethod("emfit.estep", signature(model = "map", data = "mapdata"),
  function(model, data, ufact = 1.01, eps = 1.0e-8, ...) {
    data@data$instant[is.na(data@data$counts)] <- 0
    data@data$counts[is.na(data@data$counts)] <- -1
    res <- .Call(mapfit_estep_gen_group, model, data, eps, ufact)
    list(eres=list(eb=res[[1]], ez=res[[2]], en0=res[[3]], en1=res[[4]]), llf=res[[5]])
  })

#### mstep

setMethod("emfit.mstep", signature(model = "map"),
  function(model, eres, data, stationary = TRUE, ...) {
    res <- .Call(mapfit_mstep_gen, model, eres, data)
    model@D0@x <- res[[2]]
    model@D1@x <- res[[3]]
    if (stationary)
      model@alpha <- ctmc.st(as.matrix(model@D0 + model@D1))
    else
      model@alpha <- res[[1]]
    model
  })

