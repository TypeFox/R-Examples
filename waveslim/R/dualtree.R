dualfilt1 <- function() {

  af1 <- c(0.03516384000000, 0,
           0, 0,
           -0.08832942000000, -0.11430184000000,
           0.23389032000000, 0,
           0.76027237000000, 0.58751830000000,
           0.58751830000000, -0.76027237000000,
           0, 0.23389032000000,
           -0.11430184000000, 0.08832942000000,
           0, 0,
           0, -0.03516384000000)
  af1 <- matrix(af1, ncol=2, byrow=TRUE)
  af2 <- c(0, -0.03516384000000,
           0, 0,
           -0.11430184000000, 0.08832942000000,
           0, 0.23389032000000,
           0.58751830000000, -0.76027237000000,
           0.76027237000000, 0.58751830000000,
           0.23389032000000, 0,
           -0.08832942000000, -0.11430184000000,
           0, 0,
           0.03516384000000, 0)
  af2 <- matrix(af2, ncol=2, byrow=TRUE)
  sf1 <- af1[nrow(af1):1, ]
  sf2 <- af2[nrow(af2):1, ]
  list(af = list(af1, af2), sf = list(sf1, sf2))
}

FSfarras <- function() {

  af1 <- c(0, 0,
           -0.08838834764832,  -0.01122679215254,
           0.08838834764832,   0.01122679215254,
           0.69587998903400,   0.08838834764832,
           0.69587998903400,   0.08838834764832,
           0.08838834764832,  -0.69587998903400,
           -0.08838834764832,   0.69587998903400,
           0.01122679215254,  -0.08838834764832,
           0.01122679215254,  -0.08838834764832,
           0, 0)
  af1 <- matrix(af1, ncol=2, byrow=TRUE)
  sf1 <- af1[nrow(af1):1, ]
  af2 <- c(0.01122679215254, 0,
           0.01122679215254, 0,
           -0.08838834764832, -0.08838834764832,
           0.08838834764832, -0.08838834764832,
           0.69587998903400, 0.69587998903400,
           0.69587998903400, -0.69587998903400,
           0.08838834764832, 0.08838834764832,
           -0.08838834764832, 0.08838834764832,
           0, 0.01122679215254,
           0, -0.01122679215254)
  af2 <- matrix(af2, ncol=2, byrow=TRUE)
  sf2 <- af2[nrow(af2):1, ]
  list(af = list(af1, af2), sf = list(sf1, sf2))
}

farras <- function() {

  af <- c(0, -0.01122679215254, 0, 0.01122679215254,
          -0.08838834764832, 0.08838834764832,
          0.08838834764832, 0.08838834764832,
          0.69587998903400, -0.69587998903400,
          0.69587998903400, 0.69587998903400,
          0.08838834764832, -0.08838834764832,
          -0.08838834764832, -0.08838834764832,
          0.01122679215254, 0, 0.01122679215254, 0)
  af <- matrix(af, nrow=10, byrow=TRUE)
  sf <- af[nrow(af):1, ]
  list(af = af, sf = sf)
}

cshift <- function(x, m) {

  N <- length(x)
  n <- 0:(N-1)
  n <- (n-m) %% N
  y <- x[n+1]
  y
}

afb <- function(x, af) {

  N <- length(x)
  L <- nrow(af)/2
  x <- cshift(x,-L)
  
  ## lowpass filter
  lo <- convolve(x, af[,1], conj=FALSE, type="open")
  lo <- cshift(lo,-(2*L-1))
  lo <- lo[seq(1, length(lo), by=2)]
  lo[1:L] <- lo[N/2+(1:L)] + lo[1:L]
  lo <- lo[1:(N/2)]
  
  ## highpass filter
  hi <- convolve(x, af[,2], conj=FALSE, type="open")
  hi <- cshift(hi,-(2*L-1))
  hi <- hi[seq(1, length(hi), by=2)]
  hi[1:L] <- hi[N/2+(1:L)] + hi[1:L]
  hi <- hi[1:(N/2)]

  list(lo = lo, hi = hi)
}

dualtree <- function(x, J, Faf, af) {

  ## normalization
  x <- x/sqrt(2)
  w <- vector("list", J+1)

  ## Tree 1
  w[[1]] <- vector("list", 2)
  temp <- afb(x, Faf[[1]])
  x1 <- temp$lo
  w[[1]][[1]] <- temp$hi
  if(J > 1) {
    for(j in 2:J) {
      w[[j]] <- vector("list", 2)
      temp <- afb(x1, af[[1]])
      x1 <- temp$lo
      w[[j]][[1]] <- temp$hi
    }
  }
  w[[J+1]] <- vector("list", 2)
  w[[J+1]][[1]] <- x1

  ## Tree 2
  temp <- afb(x, Faf[[2]])
  x2 <- temp$lo
  w[[1]][[2]] <- temp$hi
  if(J > 1) {
    for(j in 2:J) {
      temp <- afb(x2, af[[2]])
      x2 <- temp$lo
      w[[j]][[2]] <- temp$hi
    }
  }
  w[[J+1]][[2]] <- x2
  w
}

sfb <- function(lo, hi, sf) {

  N <- 2*length(lo)
  L <- nrow(sf)
  ## lo <- upfirdn(lo, sf[,1], 2, 1)
  lo <- c(matrix(c(rep(0, N/2), lo), nrow=2, byrow=TRUE))
  lo <- convolve(lo, sf[,1], conj=FALSE, type="open")
  lo <- cshift(lo, -L)
  ## hi <- upfirdn(hi, sf[,2], 2, 1)
  hi <- c(matrix(c(rep(0, N/2), hi), nrow=2, byrow=TRUE))
  hi <- convolve(hi, sf[,2], conj=FALSE, type="open")
  hi <- cshift(hi, -L)
  
  y <- lo + hi
  y[1:(L-2)] <- y[1:(L-2)] + y[N+1:(L-2)]
  y <- y[1:N]
  ## y = cshift(y, 1-L/2);
  y <- cshift(y, 1-L/2)
  y
}

idualtree <- function(w, J, Fsf, sf) {

  ## Tree 1
  y1 <- w[[J+1]][[1]]
  if(J > 1) {
    for(j in J:2) {
      y1 <- sfb(y1, w[[j]][[1]], sf[[1]])
    }
  }
  y1 <- sfb(y1, w[[1]][[1]], Fsf[[1]])
  
  ## Tree 2
  y2 <- w[[J+1]][[2]]
  if(J > 1) {
    for(j in J:2) {
      y2 <- sfb(y2, w[[j]][[2]], sf[[2]])
    }
  }
  y2 <- sfb(y2, w[[1]][[2]], Fsf[[2]])

  ## normalization
  y <- (y1 + y2)/sqrt(2)
  y
}

