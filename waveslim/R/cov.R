## Wavelet Variance

wave.variance <- function(x, type="eta3", p=0.025) {
  ci.gaussian <- function(x, y, p) {
    find.first <- function(v) {
      na.length <- sum(is.na(v))
      v[na.length + 1]
    }
    x.acf <- lapply(x, FUN = my.acf)
    Aj <- unlist(lapply(x.acf, FUN = function(v) sum(v*v, na.rm = TRUE))) -
      unlist(lapply(x.acf, FUN = find.first))^2 / 2
    wv.var <- 2 * Aj / unlist(lapply(x, FUN = function(v) sum(!is.na(v))))
    return(data.frame(wavevar = y, lower = y - qnorm(1-p) * sqrt(wv.var),
                      upper = y + qnorm(1-p) * sqrt(wv.var)))
  }
  
  ci.eta1 <- function(x, y, p) {
    ## x4 <- lapply(x, FUN = function(v) sum(v^4, na.rm = TRUE))
    ## eta1 <- x4.ss * unlist(lapply(x, FUN = function(v) sum(!is.na(v))))
    return(0)
  }
  
  ci.eta2 <- function(x, y, p) {
    return(0)
  }
  
  ci.eta3 <- function(x, y, p) {
    x.length <- unlist(lapply(x, FUN=function(v)sum(!is.na(v))))
    eta3 <- pmax(x.length / 2^(1:length(x)), 1)
    return(data.frame(wavevar = y, lower = eta3 * y / qchisq(1-p, eta3),
                      upper = eta3 * y / qchisq(p, eta3)))
  }

  ci.nongaussian <- function(x, y, p) {
    K <- 5
    J <- length(x)
    x.length <- unlist(lapply(x, FUN=function(v)sum(!is.na(v))))
    x.ss <- unlist(lapply(x, FUN=function(v)v[!is.na(v)]^2))
    mt.var <- numeric(J)
    for(j in 1:J) { 
      x.dpss <- dpss.taper(x.length[j], K, 4)
      V <- apply(x.dpss, 2, sum)
      J <- apply(x.dpss * x.ss[[j]], 2, sum)
      mt.var[j] <- sum((J - y[j] * V)^2) / K / x.length[j]
    }
    return(data.frame(wavevar = y, lower = y - qnorm(1-p) * sqrt(mt.var),
                      upper = y + qnorm(1-p) * sqrt(mt.var)))
  }

  x.ss <- unlist(lapply(x, FUN = function(v) sum(v*v, na.rm=TRUE)))
  x.length <- unlist(lapply(x, FUN = function(v) sum(!is.na(v))))
  y <- x.ss / x.length

  switch(type,
    gaussian = ci.gaussian(x, y, p),
    eta1 = ci.eta1(x, y, p),
    eta2 = ci.eta2(x, y, p),
    eta3 = ci.eta3(x, y, p),
    nongaussian = ci.nongaussian(x, y, p),
    stop("Invalid selection of \"type\" for the confidence interval"))
}

##plot.var <- function(x, y=NA, ylim=range(x, y, na.rm=TRUE))
##{
##  n <- dim(x)[1]
##  plot(2^(0:(n-1)), x[,1], axes=FALSE, type="n", log="xy", ylim=ylim)
##  axis(1, at=2^(0:(n-1)))
##  axis(2)
##  polyci(x[,1], x[,-1], -1)
##  if(any(!is.na(y))) { polyci(y[,1], y[,-1], 1, color=5) }
##  abline(h=0, lty=2)
##}

## Wavelet Covariance

wave.covariance <- function(x, y) {
  my.acf.na <- function(v) {
    v <- v[!is.na(v)]
    my.acf(v)
  }
  my.ccf.na <- function(u, v) {
    u <- u[!is.na(u)]
    v <- v[!is.na(v)]
    n <- length(u)
    u <- c(u, rep(0, n))
    v <- c(v, rep(0, n))
    n <- length(u)
    x <- Re(fft(fft(u) * Conj(fft(v)), inverse=TRUE)) / 2 / n^2
    x[c((n %/% 2):n, 1:(n %/% 2 - 1))]
  }
  compute.sum.xy.ccvs <- function(x, y) {
    l <- length(x)
    xy <- numeric(l)
    for(i in 1:l)
      xy[i] <- sum(my.ccf.na(x[[i]], y[[i]])^2)
    xy
  }
  compute.xy.acvs <- function(x, y) {
    l <- length(x)
    xy <- vector("list", l)
    for(i in 1:l) {
      z <- x[[i]] * y[[i]]
      xy[[i]] <- c(rev(z), z[-1])
    }
    xy
  }
  per <- function (z) {
    n <- length(z)
    (Mod(fft(z))^2/n)[1:(n%/%2 + 1)]
  }
  per2 <- function(x, y) {
    n <- length(x)
    fft.x <- fft(x)
    fft.y <- fft(y)
    ((Conj(fft.x) * fft.y)/n)[1:(n %/% 2 + 1)]
  }
  
  l <- length(x)
  xy <- vector("list", l)
  for(i in 1:l)
    xy[[i]] <- as.vector(x[[i]] * y[[i]])
  z.ss <- unlist(lapply(xy, sum, na.rm=TRUE))
  x.na <- lapply(x, is.na)
  for(i in 1:l)
    x.na[[i]] <- !x.na[[i]]
  z.length <- unlist(lapply(x.na, sum))

  zz <- z.ss / z.length
  names(zz) <- names(x)

  x.acvs <- lapply(x, my.acf.na)
  y.acvs <- lapply(y, my.acf.na)
  sum.xy.acvs <- unlist(lapply(compute.xy.acvs(x.acvs, y.acvs), sum))
  sum.squared.xy.ccvs <- compute.sum.xy.ccvs(x, y)
  var.gamma <- (sum.xy.acvs + sum.squared.xy.ccvs) / 2 / z.length

  out <- data.frame(wavecov = zz, lower = zz - qnorm(.975) * sqrt(var.gamma),
                    upper = zz + qnorm(.975) * sqrt(var.gamma))
  return(as.matrix(out))
}

##polyci <- function(x, xci, sp, color=2)
##{
##  n <- length(x)
##  y <- 2^(0:(n-1)+sp*.05)
##  delta <- y - 2^(0:(n-1))
##  for(i in 1:n){
##     polygon(c(y[i] + .6*delta[i], y[i] + .6*delta[i], y[i] - .6*delta[i],
##        y[i] - .6*delta[i]), c(xci[i,], xci[i,2:1]), border=FALSE,
##        col=color, lty=1)
##  }
##  points(y, x, pch="-")
##}

##plot.cov <- function(x, ylim=range(x,0))
##{
##  n <- dim(x)[1]
##  plot(2^(0:(n-1)), x[,1], axes=FALSE, type="n", log="x", ylim=ylim)
##  axis(1, at=2^(0:(n-1)))
##  axis(2)
##  polyci(x[,1], x[,-1], 1)
##  abline(h=0, lty=2)
##}

## Wavelet Correlation

wave.correlation <- function(x, y, N, p = .975) {
  sum.of.squares <- function(x) { sum(x^2, na.rm=TRUE) / sum(!is.na(x)) }
  sum.of.not.squares <- function(x) { sum(x, na.rm=TRUE) / sum(!is.na(x)) }

  l <- length(x)
  xy <- vector("list", l); xy.abs <- vector("list", l)
  for(i in 1:l) {
    xy[[i]] <- as.vector(x[[i]] * y[[i]])
    xy.abs[[i]] <- as.vector(abs(x[[i]] * y[[i]]))
  }
  xy.cov <- unlist(lapply(xy, sum.of.not.squares))
  
  x.var <- unlist(lapply(x, sum.of.squares))
  y.var <- unlist(lapply(y, sum.of.squares))
  
  xy.cor <- xy.cov / sqrt(x.var * y.var)
  n <- trunc(N/2^(1:l))
  out <- data.frame(wavecor=xy.cor,
                    lower=tanh(atanh(xy.cor)-qnorm(p)/sqrt(n-3)),
                    upper=tanh(atanh(xy.cor)+qnorm(p)/sqrt(n-3)))
  return(as.matrix(out))
}

##plot.cor <- function(x, ylim=c(-1,1), cex=NULL)
##{
##  n <- dim(x)[1]
##  plot(2^(0:(n-1)), x[,1], axes=FALSE, type="n", log="x", ylim=ylim, cex=cex)
##  axis(1, at=2^(0:(n-1)), cex=cex)
##  axis(2, cex=cex)
##  polyci(x[,1], x[,-1], 1)
##  abline(h=0, lty=2)
##}

## Plotting functions for wavelet variance and covariance

## Wavelet cross-covariance 

spin.covariance <- function(x, y, lag.max = NA) {
  xx <- zz <- x[!is.na(x)]
  yy <- y[!is.na(y)]
  n.length <- length(xx)
  xx.length <- min(length(xx)-1, lag.max, na.rm=TRUE)

  lag1 <- numeric(xx.length + 1)
  lag2 <- numeric(xx.length + 1)
  for(i in 1:(xx.length+1)) {
      lag1[i] <- sum(xx * yy, na.rm=TRUE) / n.length
      lag2[i] <- sum(zz * yy, na.rm=TRUE) / n.length
      xx <- c(xx[2:n.length], NA)
      zz <- c(NA, zz[1:(n.length-1)])
    }
  c(rev(lag2[-1]), lag1)
}

spin.correlation <- function(x, y, lag.max = NA) {
  xx <- zz <- x[!is.na(x)]
  yy <- y[!is.na(y)]
  n.length <- length(xx)
  xx.length <- min(length(xx)-1, lag.max, na.rm=TRUE)
  xx.var <- mean(xx^2)
  yy.var <- mean(yy^2)

  lag1 <- numeric(xx.length + 1)
  lag2 <- numeric(xx.length + 1)
  for(i in 1:(xx.length+1)) {
      lag1[i] <- sum(xx * yy, na.rm=TRUE) / sqrt(xx.var * yy.var) / n.length
      lag2[i] <- sum(zz * yy, na.rm=TRUE) / sqrt(xx.var * yy.var) / n.length
      xx <- c(xx[2:n.length], NA)
      zz <- c(NA, zz[1:(n.length-1)])
    }
  c(rev(lag2[-1]), lag1)
}

##edof <- function(x) {
##  x <- x[!is.na(x)]
##  n <- length(x)
##  x.acf <- my.acf(x)
##  n * x.acf[1]^2 /
##    sum((1 - abs(seq(-n+1,n-1))/n) * c(rev(x.acf[-1]), x.acf)^2)
##}
