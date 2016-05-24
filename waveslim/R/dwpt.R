dwpt <- function(x, wf="la8", n.levels=4, boundary="periodic") {
  N <- length(x)
  J <- n.levels
  if(N/2^J != trunc(N/2^J))
    stop("Sample size is not a power of 2")
  if(2^J > N)
    stop("wavelet transform exceeds sample size in dwt")

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  y <- vector("list", sum(2^(1:J)))
  crystals1 <- rep(1:J, 2^(1:J))
  crystals2 <- unlist(apply(as.matrix(2^(1:J) - 1), 1, seq, from=0))
  names(y) <- paste("w", crystals1, ".", crystals2, sep="")

  for(j in 1:J) {
    jj <- min((1:length(crystals1))[crystals1 == j])
    for(n in 0:(2^j/2-1)) {
      if(j > 1)
        x <- y[[(1:length(crystals1))[crystals1 == j-1][n+1]]]
      W <- V <- numeric(N/2^j)
      if(n %% 2 == 0) {
        z <- .C("dwt", as.double(x), as.integer(N/2^(j-1)), L, h, g, 
	  W=as.double(W), V=as.double(V), PACKAGE="waveslim")
        y[[jj + 2*n + 1]] <- z$W
        y[[jj + 2*n]] <- z$V
      }
      else {
        z <- .C("dwt", as.double(x), as.integer(N/2^(j-1)), L, h, g,
                W=as.double(W), V=as.double(V), PACKAGE="waveslim")
        y[[jj + 2*n]] <- z$W
        y[[jj + 2*n + 1 ]] <- z$V
      }
    }
  }
  attr(y, "wavelet") <- wf
  return(y)
}

idwpt <- function(y, y.basis)
{
  J <- trunc(log(length(y), 2))

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  for(j in J:1) {
    a <- min((1:length(rep(1:J, 2^(1:J))))[rep(1:J, 2^(1:J)) == j])
    b <- max((1:length(rep(1:J, 2^(1:J))))[rep(1:J, 2^(1:J)) == j])
    n <- a
    while(n <= b) {
      if(y.basis[n]) {
        m <- length(y[[n]])
        XX <- numeric(2 * m)
        if(floor((n-a)/2) %% 2 == 0)
          X <- .C("idwt", as.double(y[[n+1]]), as.double(y[[n]]),
                  as.integer(m), L, h, g, out=as.double(XX),
                  PACKAGE="waveslim")$out
        else
          X <- .C("idwt", as.double(y[[n]]), as.double(y[[n+1]]), 
                  as.integer(m), L, h, g, out=as.double(XX),
                  PACKAGE="waveslim")$out
        if(j != 1) {
          y[[a-(b-a+1)/2 + (n-a)/2]] <- X
          y.basis[[a-(b-a+1)/2 + (n-a)/2]] <- 1
        }
        n <- n + 2
      }
      else { n <- n + 1 }
    }
  }
  return(X)
}

##plot.dwpt <- function(x, n.levels, pgrid=TRUE)
##{
##  J <- n.levels
##  scales <- rep(1:J, 2^(1:J))
##  y <- matrix(NA, 2*length(x[[1]]), J)
##  for(j in 1:J) {
##    a <- min((1:length(scales))[scales == j])
##    b <- max((1:length(scales))[scales == j])
##    y[, j] <- unlist(x[a:b])
##    x.length <- length(y[, j])
##  }
##  plot(ts(y), ylim=c(-.45,.45))
##  if(pgrid) {
##    lines(x.length * c(0,1), c(0,0), lty=2)
##    for(j in 1:J) {
##      lines(x.length * c(0,1), c(-j,-j), lty=2)
##      for(n in 0:2^j) lines(x.length * c(n/2^j, n/2^j), c(-j,-(j-1)), lty=2)
##    }
##  }
##  title(ylab="Level")
##}

basis <- function(x, basis.names)
{
  m <- length(x)
  n <- length(basis.names)
  y <- numeric(m)
  for(i in 1:n) { y <- y + as.integer(names(x) == basis.names[i]) }
  return(y)
}

ortho.basis <- function(xtree) {
  J <- trunc(log(length(xtree), 2))
  X <- vector("list", J)
  X[[1]] <- xtree[rep(1:J, 2^(1:J)) == 1]
  for(i in 2:J) {
    for(j in i:J) {
      if(i == 2) X[[j]] <- xtree[rep(1:J, 2^(1:J)) == j]
        X[[j]] <- X[[j]] + 2 * c(apply(matrix(xtree[rep(1:J, 2^(1:J)) == i-1]),
                                       1, rep, 2^(j-i+1)))
    }
  }
  X[[J]][X[[J]] == 0] <- 1
  ifelse(unlist(X) == 1, 1, 0)
}

##plot.basis <- function(xtree)
##{
##  J <- trunc(log(length(xtree), base=2))
##  j <- rep(1:J, 2^(1:J))
##  n <- unlist(apply(matrix(2^(1:J)-1), 1, seq, from=0))
##  basis <- ifelse(xtree, paste("w", j, ".", n, sep=""), NA)
##  pgrid.plot(basis[basis != "NA"])
##  invisible()
##}

phase.shift.packet <- function(z, wf, inv=FALSE)
{
  ## Center of energy
  coe <- function(g)
    sum(0:(length(g)-1) * g^2) / sum(g^2)

  J <- length(x) - 1
  g <- wave.filter(wf)$lpf
  h <- wave.filter(wf)$hpf

  if(!inv) {
    for(j in 1:J) {
      ph <- round(2^(j-1) * (coe(g) + coe(h)) - coe(g), 0)
      Nj <- length(x[[j]])
      x[[j]] <- c(x[[j]][(ph+1):Nj], x[[j]][1:ph])
    }
    ph <- round((2^J-1) * coe(g), 0)
    J <- J + 1
    x[[J]] <- c(x[[J]][(ph+1):Nj], x[[J]][1:ph])
  } else {
    for(j in 1:J) {
      ph <- round(2^(j-1) * (coe(g) + coe(h)) - coe(g), 0)
      Nj <- length(x[[j]])
      x[[j]] <- c(x[[j]][(Nj-ph+1):Nj], x[[j]][1:(Nj-ph)])
    }
    ph <- round((2^J-1) * coe(g), 0)
    J <- J + 1
    x[[J]] <- c(x[[j]][(Nj-ph+1):Nj], x[[j]][1:(Nj-ph)])
  }
  return(x)
}

modwpt <- function(x, wf="la8", n.levels=4, boundary="periodic")
{
  N <- length(x); storage.mode(N) <- "integer"
  J <- n.levels
  if(2^J > N) stop("wavelet transform exceeds sample size in modwt")

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  ht <- dict$hpf/sqrt(2)
  storage.mode(ht) <- "double"
  gt <- dict$lpf/sqrt(2)
  storage.mode(gt) <- "double"

  y <- vector("list", sum(2^(1:J)))
  yn <- length(y)
  crystals1 <- rep(1:J, 2^(1:J))
  crystals2 <- unlist(apply(as.matrix(2^(1:J) - 1), 1, seq, from=0))
  names(y) <- paste("w", crystals1, ".", crystals2, sep="")

  W <- V <- numeric(N)
  storage.mode(W) <- storage.mode(V) <- "double"
  for(j in 1:J) {
    index <- 0
    jj <- min((1:yn)[crystals1 == j])
    for(n in 0:(2^j / 2 - 1)) {
      index <- index + 1
      if(j > 1)
        x <- y[[(1:yn)[crystals1 == j-1][index]]]
      if(n %% 2 == 0) {
        z <- .C("modwt", as.double(x), N, as.integer(j), L, ht, gt, 
                W = W, V = V, PACKAGE = "waveslim")[7:8]
        y[[jj + 2*n + 1]] <- z$W
        y[[jj + 2*n]] <- z$V
      }
      else {
        z <- .C("modwt", as.double(x), N, as.integer(j), L, ht, gt, 
                W = W, V = V, PACKAGE = "waveslim")[7:8]
        y[[jj + 2*n]] <- z$W
        y[[jj + 2*n + 1 ]] <- z$V
      }
    }
  }
  attr(y, "wavelet") <- wf
  return(y)
}

dwpt.brick.wall <- function(x, wf, n.levels, method="modwpt")
{
  N <- length(x[[1]])
  m <- wave.filter(wf)$length
  J <- n.levels
  crystals1 <- rep(1:J, 2^(1:J))
  crystals2 <- unlist(apply(as.matrix(2^(1:J) - 1), 1, seq, from=0))

  if(method=="dwpt") {
    ## for DWPT
    for(j in 1:J) {
      jj <- min((1:length(crystals1))[crystals1 == j])
      L <- switch(j,
                  (m-2)/2,
                  ((m-2)/2 + floor(m/4)),
                  ((m-2)/2 + floor((m/2 + floor(m/4))/2)))
      if(is.null(L)) L <- (m-2)
      for(n in 0:(2^j-1))
        x[[jj+n]][1:L] <- NA
    }
  }
  else {
    ## for MODWPT
    for(j in 1:J) {
      jj <- min((1:length(crystals1))[crystals1 == j])
      L <- min((2^j - 1) * (m - 1), N)
      for(n in 0:(2^j-1))
        x[[jj+n]][1:L] <- NA
    }
  }
  return(x)
}

css.test <- function(y) 
{
  K <- length(y)
  test <- numeric(K)

  for(k in 1:K) {
    x <- y[[k]]
    x <- x[!is.na(x)]
    n <- length(x)
    plus <- 1:n/(n - 1) - cumsum(x^2)/sum(x^2)
    minus <- cumsum(x^2)/sum(x^2) - 0:(n - 1)/(n - 1)
    D <- max(abs(plus), abs(minus))
    if(D < 1.224/(sqrt(n) + 0.12 + 0.11/sqrt(n))) test[k] <- 1
  }
  return(test)
}

entropy.test <- function(y)
{
  K <- length(y)
  test <- numeric(K)

  for(k in 1:K) {
    x <- y[[k]]
    test[k] <- sum(x^2 * log(x^2), na.rm=TRUE)
  }
  return(test)
}

cpgram.test <- function(y, p=0.05, taper=0.1)
{
  K <- length(y)
  test <- numeric(K)
  
  for(k in 1:K) {
    x <- y[[k]]
    x <- x[!is.na(x)]
    x <- spec.taper(scale(x, center=TRUE, scale=FALSE), p=taper)
    y <- Mod(fft(x))^2/length(x)
    y[1] <- 0
    n <- length(x)
    x <- (0:(n/2))/n
    if(length(x) %% 2 == 0) {
      n <- length(x) - 1
      y <- y[1:n]
      x <- x[1:n]
    }
    else y <- y[1:length(x)]
    mp <- length(x) - 1
    if(p == 0.05)
      crit <- 1.358/(sqrt(mp) + 0.12 + 0.11/sqrt(mp))
    else {
      if(p == 0.01) crit <- 1.628/(sqrt(mp) + 0.12 + 0.11/sqrt(mp))
      else stop("critical value is not known")
    }
    D <- abs(cumsum(y)/sum(y) -  0:mp/mp)
    if(max(D) < crit) test[k] <- 1
  }
  return(test)
}

portmanteau.test <- function(y, p = 0.05, type = "Box-Pierce")
{
  K <- length(y)
  test <- numeric(K)

  for(k in 1:K) {
  x <- y[[k]]
  x <- x[!is.na(x)]
  n <- length(x)
  h <- trunc(n/2)
  x.acf <- my.acf(x)[1:(h+1)]
  x.acf <- x.acf / x.acf[1];
  if(type == "Box-Pierce")
    test[k] <- ifelse(n * sum((x.acf[-1])^2) > qchisq(1-p, h), 0, 1)
  else 
    test[k] <- ifelse(n*(n+2) * sum((x.acf[-1])^2 / (n - h:1)) > 
                      qchisq(1-p, h), 0, 1)
  }
  return(test)
}

