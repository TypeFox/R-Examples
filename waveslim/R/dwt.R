dwt <- function(x, wf="la8", n.levels=4, boundary="periodic")
{
  switch(boundary,
    "reflection" =  x <- c(x, rev(x)),
    "periodic" = invisible(),
    stop("Invalid boundary rule in dwt"))
  N <- length(x)
  J <- n.levels
  if(N/2^J != trunc(N/2^J))
    stop("Sample size is not divisible by 2^J")
  if(2^J > N)
    stop("wavelet transform exceeds sample size in dwt")

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  y <- vector("list", J+1)
  names(y) <- c(paste("d", 1:J, sep=""), paste("s", J, sep=""))
  for(j in 1:J) {
    W <- V <- numeric(N/2^j)
    out <- .C("dwt", as.double(x), as.integer(N/2^(j-1)), L, h, g, 
              W=as.double(W), V=as.double(V), PACKAGE="waveslim")[6:7]
    y[[j]] <- out$W
    x <- out$V
  }
  y[[J+1]] <- x
  class(y) <- "dwt"
  attr(y, "wavelet") <- wf
  attr(y, "boundary") <- boundary
  return(y)
}

dwt.nondyadic <- function(x)
{
  M <- length(x)
  N <- 2^(ceiling(log(M, 2)))
  xx <- c(x, rep(0, N - M))
  y <- dwt(xx)
  
  J <- length(y) - 1
  for(j in 1:J)
    y[[j]] <- y[[j]][1:trunc(M/2^j)]
  return(y)
}

idwt <- function(y)
{
  ctmp <- class(y)
  if(is.null(ctmp) || all(ctmp != "dwt"))
    stop("argument `y' is not of class \"dwt\"")

  J <- length(y) - 1

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  jj <- paste("s", J, sep="")
  X <- y[[jj]]
  for(j in J:1) {
    jj <- paste("d", j, sep="")
    N <- length(X)
    XX <- numeric(2 * length(y[[jj]]))
    X <- .C("idwt", as.double(y[[jj]]), as.double(X), as.integer(N), L, 
            h, g, out=as.double(XX), PACKAGE="waveslim")$out
  }
  if(attr(y, "boundary") == "reflection") return(X[1:N])
  else return(X)
}

modwt <- function(x, wf="la8", n.levels=4, boundary="periodic")
{
  switch(boundary,
    "reflection" =  x <- c(x, rev(x)),
    "periodic" = invisible(),
    stop("Invalid boundary rule in modwt"))
  N <- length(x)
  storage.mode(N) <- "integer"
  J <- n.levels
  if(2^J > N)
    stop("wavelet transform exceeds sample size in modwt")

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  ht <- dict$hpf / sqrt(2)
  storage.mode(ht) <- "double"
  gt <- dict$lpf / sqrt(2)
  storage.mode(gt) <- "double"

  y <- vector("list", J+1)
  names(y) <- c(paste("d", 1:J, sep=""), paste("s", J, sep=""))
  W <- V <- numeric(N)
  storage.mode(W) <- "double"
  storage.mode(V) <- "double"
  
  for(j in 1:J) {
    out <- .C("modwt", as.double(x), N, as.integer(j), L, ht, gt, 
              W=W, V=V, PACKAGE="waveslim")[7:8]
    y[[j]] <- out$W
    x <- out$V
  }
  y[[J+1]] <- x
  class(y) <- "modwt"
  attr(y, "wavelet") <- wf
  attr(y, "boundary") <- boundary
  return(y)
}

imodwt <- function(y)
{
  ctmp <- class(y)
  if(is.null(ctmp) || all(ctmp != "modwt"))
    stop("argument `y' is not of class \"modwt\"")

  J <- length(y) - 1

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  ht <- dict$hpf / sqrt(2)
  storage.mode(ht) <- "double"
  gt <- dict$lpf / sqrt(2)
  storage.mode(gt) <- "double"

  jj <- paste("s", J, sep="")
  X <- y[[jj]]
  N <- length(X)
  storage.mode(N) <- "integer"
  XX <- numeric(N)
  storage.mode(XX) <- "double"
  for(j in J:1) {
    jj <- paste("d", j, sep="")
    X <- .C("imodwt", as.double(y[[jj]]), as.double(X), N, as.integer(j), 
            L, ht, gt, out=XX, PACKAGE="waveslim")$out
  }
  if(attr(y, "boundary") == "reflection") return(X[1:(N/2)])
  else return(X)
}

brick.wall <- function(x, wf, method="modwt")
{
  m <- wave.filter(wf)$length
  for(j in 1:(length(x)-1)) {
    if(method == "dwt")
      n <- ceiling((m - 2) * (1 - 1/2^j))
    else
      n <- (2^j - 1) * (m - 1)
    n <- min(n, length(x[[j]]))
    x[[j]][1:n] <- NA
  }
  x[[j+1]][1:n] <- NA
  return(x)
}

phase.shift <- function(z, wf, inv=FALSE)
{
  coe <- function(g)
    sum(0:(length(g)-1) * g^2) / sum(g^2)

  J <- length(z) - 1
  g <- wave.filter(wf)$lpf
  h <- wave.filter(wf)$hpf

  if(!inv) {
    for(j in 1:J) {
      ph <- round(2^(j-1) * (coe(g) + coe(h)) - coe(g), 0)
      Nj <- length(z[[j]])
      z[[j]] <- c(z[[j]][(ph+1):Nj], z[[j]][1:ph])
    }
    ph <- round((2^J-1) * coe(g), 0)
    J <- J + 1
    z[[J]] <- c(z[[J]][(ph+1):Nj], z[[J]][1:ph])
  } else {
    for(j in 1:J) {
      ph <- round(2^(j-1) * (coe(g) + coe(h)) - coe(g), 0)
      Nj <- length(z[[j]])
      z[[j]] <- c(z[[j]][(Nj-ph+1):Nj], z[[j]][1:(Nj-ph)])
    }
    ph <- round((2^J-1) * coe(g), 0)
    J <- J + 1
    z[[J]] <- c(z[[J]][(Nj-ph+1):Nj], z[[J]][1:(Nj-ph)])
  }
  return(z)
}

mra <- function(x, wf="la8", J=4, method="modwt", boundary="periodic")
{
  switch(boundary,
         "reflection" =  x <- c(x, rev(x)),
         "periodic" = invisible(),
         stop("Invalid boundary rule in mra"))
  n <- length(x)
  
  if(method == "modwt")
    x.wt <- modwt(x, wf, J, "periodic")
  else
    x.wt <- dwt(x, wf, J, "periodic")
  x.mra <- vector("list", J+1)
  
  ## Smooth
  zero <- vector("list", J+1)
  names(zero) <- c(paste("d", 1:J, sep=""), paste("s", J, sep=""))
  class(zero) <- method
  attr(zero, "wavelet") <- wf
  attr(zero, "boundary") <- boundary
  zero[[J+1]] <- x.wt[[J+1]]
  if(method == "modwt") {
    for(k in 1:J)
      zero[[k]] <- numeric(n)
    x.mra[[J+1]] <- imodwt(zero)
  } else {
    for(k in 1:J)
      zero[[k]] <- numeric(n/2^k)
    x.mra[[J+1]] <- idwt(zero)
  }

  ## Details
  for(j in J:1) {
    zero <- vector("list", j+1)
    names(zero) <- c(paste("d", 1:j, sep=""), paste("s", j, sep=""))
    class(zero) <- method
    attr(zero, "wavelet") <- wf
    attr(zero, "boundary") <- boundary
    zero[[j]] <- x.wt[[j]]
    if(method == "modwt") {
      if(j != 1) {
        for(k in c(j+1,(j-1):1))
          zero[[k]] <- numeric(n)
      } else {
        zero[[j+1]] <- numeric(n)
      }
      x.mra[[j]] <- imodwt(zero)
    } else {
      zero[[j+1]] <- numeric(n/2^j)
      if(j != 1) {
        for(k in (j-1):1)
          zero[[k]] <- numeric(n/2^k)
      }
      x.mra[[j]] <- idwt(zero)
    }
  }

  names(x.mra) <- c(paste("D", 1:J, sep=""), paste("S", J, sep=""))
  if(boundary == "reflection") { 
    for(j in (J+1):1)
      x.mra[[j]] <- x.mra[[j]][1:(n/2)]
    return(x.mra)
  } else {
    return(x.mra)
  }
}

