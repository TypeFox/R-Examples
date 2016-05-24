dualtree2D <- function(x, J, Faf, af) {
  ## normalization
  x <- x/sqrt(2)
  w <- vector("list", J+1)
  ## Tree 1
  w[[1]] <- vector("list", 2)
  temp <- afb2D(x, Faf[[1]])              # stage 1
  x1 <- temp$lo
  w[[1]][[1]] <- temp$hi
  if (J > 1) {
    for (j in 2:J) {
      w[[j]] <- vector("list", 2)
      temp <- afb2D(x1, af[[1]])          # remaining stages
      x1 <- temp$lo
      w[[j]][[1]] <- temp$hi
    }
  }
  w[[J+1]] <- vector("list", 2)
  w[[J+1]][[1]] <- x1                     # lowpass subband
  ## Tree 2
  temp <- afb2D(x, Faf[[2]])              # stage 1
  x2 <- temp$lo
  w[[1]][[2]] <- temp$hi
  if (J > 1) {
    for (j in 2:J) {
      temp <- afb2D(x2, af[[2]])          # remaining stages
      x2 <- temp$lo
      w[[j]][[2]] <- temp$hi
    }
  }
  w[[J+1]][[2]] <- x2                     # lowpass subband
  ## sum and difference
  for (j in 1:J) {
    for (m in 1:3) {
      A <- w[[j]][[1]][[m]]
      B <- w[[j]][[2]][[m]]
      w[[j]][[1]][[m]] <- (A + B) / sqrt(2)
      w[[j]][[2]][[m]] <- (A - B) / sqrt(2)
    }
  }
  return(w)
}

afb2D <- function(x, af1, af2=NULL) {
  if (is.null(af2)) {
    af2 <- af1
  }
  ## filter along columns
  temp <- afb2D.A(x, af1, 1)
  L <- temp$lo
  H <- temp$hi
  ## filter along rows
  hi <- vector("list", 3)
  temp <- afb2D.A(L, af2, 2)
  lo <- temp$lo
  hi[[1]] <- temp$hi
  temp <- afb2D.A(H, af2, 2)
  hi[[2]] <- temp$lo
  hi[[3]] <- temp$hi
  list(lo = lo, hi = hi)
}

afb2D.A <- function(x, af, d) {
  lpf <- af[,1]     # lowpass filter
  hpf <- af[,2]     # highpass filter
  if (d == 2) {
    x <- t(x)
  }
  ## x <- matrix(1:32, 32, 64)
  N <- nrow(x)
  L <- nrow(af) / 2
  x <- cshift2D(x, -L)
  ## image(x, col=rainbow(16))

  ## lo <- upfirdn(x, lpf, 1, 2)
  lo <- convolve2D(x, lpf, conj=FALSE, type="open")
  lo <- cshift2D(lo, -(2 * L - 1))
  lo <- lo[seq(1, nrow(lo), by=2),]
  lo[1:L,] <- lo[1:L,] + lo[1:L + N/2,]
  lo <- lo[1:(N/2),]

  ## hi <- upfirdn(x, hpf, 1, 2)
  hi <- convolve2D(x, hpf, conj=FALSE, type="open")
  hi <- cshift2D(hi, -(2 * L - 1))
  hi <- hi[seq(1, nrow(hi), by=2),]
  hi[1:L,] <- hi[1:L,] + hi[1:L + N/2,]
  hi <- hi[1:(N/2),]

  if (d == 2) {
    lo <- t(lo)
    hi <- t(hi)
  }
  list(lo = lo, hi = hi)
}

cshift2D <- function(x, m) {
  N <- nrow(x)
  n <- 0:(N-1)
  n <- (n-m) %% N
  y <- x[n+1,]
  return(y)
}

convolve2D <- function(x, y, conj=TRUE, type=c("circular", "open")) {
  ## Generalize convolve to handle vector arrays by calling mvfft()
  type <- match.arg(type)
  n <- nrow(x)
  ny <- length(y)
  Real <- is.numeric(x) && is.numeric(y)
  if (type == "circular") {
    if (ny != n) {
      stop("length mismatch in convolution")
    }
  } else {
    n1 <- ny - 1
    x <- rbind(matrix(0, n1, ncol(x)), x)
    n <- length(y <- c(y, rep.int(0, n - 1)))
  }
  x <- mvfft(mvfft(x) * (if (conj) Conj(fft(y)) else fft(y)), inverse=TRUE)
  (if (Real) Re(x) else x) / n
}

idualtree2D <- function(w, J, Fsf, sf) {
  ## sum and difference
  for (k in 1:J) {
    for (m in 1:3) {
      A <- w[[k]][[1]][[m]]
      B <- w[[k]][[2]][[m]]
      w[[k]][[1]][[m]] <- (A+B)/sqrt(2)
      w[[k]][[2]][[m]] <- (A-B)/sqrt(2)
    }
  }
  ## Tree 1
  y1 <- w[[J+1]][[1]]
  if (J > 1) {
    for (j in J:2) {
      y1 <- sfb2D(y1, w[[j]][[1]], sf[[1]])
    }
  }
  y1 <- sfb2D(y1, w[[1]][[1]], Fsf[[1]])
  ## Tree 2
  y2 <- w[[J+1]][[2]]
  if (J > 1) {
    for (j in J:2) {
      y2 <- sfb2D(y2, w[[j]][[2]], sf[[2]])
    }
    y2 <- sfb2D(y2, w[[1]][[2]], Fsf[[2]])
  }
  ## normalization
  y <- (y1 + y2)/sqrt(2)
  return(y)
}

sfb2D <- function(lo, hi, sf1, sf2=NULL) {
  if (is.null(sf2)) {
    sf2 <- sf1
  }
  ## filter along rows
  lo <- sfb2D.A(lo, hi[[1]], sf2, 2)
  hi <- sfb2D.A(hi[[2]], hi[[3]], sf2, 2)
  ## filter along columns
  y <- sfb2D.A(lo, hi, sf1, 1)
  return(y)
}

sfb2D.A <- function(lo, hi, sf, d) {
  lpf <- sf[,1]     # lowpass filter
  hpf <- sf[,2]     # highpass filter
  if (d == 2) {
    lo <- t(lo)
    hi <- t(hi)
  }
  N <- 2 * nrow(lo)
  M <- ncol(lo)
  L <- nrow(sf)
  ## y = upfirdn(lo, lpf, 2, 1) + upfirdn(hi, hpf, 2, 1);
  lo <- c(matrix(c(rep(0, length(lo)), c(lo)), nrow=2, byrow=TRUE))
  lo <- matrix(lo, N, M)
  lo <- convolve2D(lo, lpf, conj=FALSE, type="open")
  lo <- cshift2D(lo, -L)
  hi <- c(matrix(c(rep(0, length(hi)), c(hi)), nrow=2, byrow=TRUE))
  hi <- matrix(hi, N, M)
  hi <- convolve2D(hi, hpf, conj=FALSE, type="open")
  hi <- cshift2D(hi, -L)
  y <- lo + hi
  y[1:(L-2),] <- y[1:(L-2),] + y[N+1:(L-2),]
  y <- y[1:N,]
  y <- cshift2D(y, 1 - L/2)
  if (d == 2) {
    y <- t(y)
  }
  return(y)
}

pm <- function(a, b) {
  u <- (a + b) / sqrt(2)
  v <- (a - b) / sqrt(2)
  list(u=u, v=v)
}

