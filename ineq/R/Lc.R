Lc <- function(x, n = rep(1, length(x)), plot = FALSE)
{
    ina <- !is.na(x)    
    n <- n[ina]
    x <- as.numeric(x)[ina]
    k <- length(x)
    o <- order(x)
    x <- x[o]
    n <- n[o]
    x <- n*x
    p <- cumsum(n)/sum(n)
    L <- cumsum(x)/sum(x)
    p <- c(0,p)
    L <- c(0,L)
    L2 <- L * mean(x)/mean(n)
    Lc <- list(p,L,L2)
    names(Lc) <- c("p", "L", "L.general")
    class(Lc) <- "Lc"
    if(plot) plot(Lc)
    Lc
}

plot.Lc <- function(x, general=FALSE, lwd=2,xlab="p",ylab="L(p)",
  main="Lorenz curve", las=1, ...)
{
    if(!general)
      L <- x$L
    else
      L <- x$L.general
    plot(x$p, L, type="l", main=main, lwd=lwd, xlab=xlab, ylab=ylab, xaxs="i",
      yaxs="i", las=las, ...)
    abline(0,max(L))
}

lines.Lc <- function(x, general=FALSE, lwd=2, ...)
{
    if(!general)
      L <- x$L
    else
      L <- x$L.general
    lines(x$p, L, lwd=lwd, ...)
}

plot.theorLc <- function(x, parameter=NULL, xlab="p", ylab="L(p)", lwd=2, las=1, ...)
{
  dummy <- (0:1000)*0.001
  if(is.null(parameter))
    plot(dummy, x(dummy), type="l", xlab=xlab, ylab=ylab, xaxs="i", yaxs="i",
      lwd=lwd, las=las, ...)
  else
    plot(dummy, x(dummy,parameter=parameter), type="l", xlab=xlab,
      ylab=ylab, xaxs="i", yaxs="i", lwd=lwd, las=las, ...)
  abline(0,1)
}

lines.theorLc <- function(x, parameter=NULL, lwd=2, col=2, ...)
{
  dummy <- (0:1000)*0.001
  if(is.null(parameter))
    lines(dummy, x(dummy), lwd=lwd, col=col, ...)
  else
    lines(dummy, x(dummy,parameter=parameter), lwd=lwd, col=col, ...)
}

theorLc <- function(type=c("Singh-Maddala", "Dagum", "lognorm", "Pareto",
    "exponential"), parameter=0)
{
  switch(match.arg(type),
  "Singh-Maddala" = rval <- function(p) {Lc.singh(p,parameter=parameter)},
  Dagum = rval <- function(p) {Lc.dagum(p,parameter=parameter)},
  lognorm = rval <- function(p) {Lc.lognorm(p,parameter=parameter)},
  Pareto = rval <- function(p) {Lc.pareto(p,parameter=parameter)},
  exponential = rval <- function(p) {Lc.exp(p)})
  class(rval) <- "theorLc"
  return(rval)
}

Lc.exp <- function(p)
{
  elc <- 1/(1-p)
  elc <- (1-p)*log(elc)
  elc <- p - elc
  elc
}
class(Lc.exp) <- "theorLc"

Lc.lognorm <- function(p, parameter=1)
{
  if(parameter[1]>0)
    sigma <- parameter[1]
  else
  {
    warning("inadmissible parameter. default parameter=1 is used.")
    sigma <- 1
  }
  loglc <- p
  loglc[!loglc==0 & !loglc==1] <- pnorm(qnorm(loglc[!loglc==0 & !loglc==1]) - sigma)
  loglc
}
class(Lc.lognorm) <- "theorLc"

Lc.pareto <- function(p, parameter=2)
{
  if(parameter[1]>1) k<-(parameter[1]-1)/parameter[1]
  else
  {
    warning("inadmissible parameter. default parameter=2 is used.")
    k <- 0.5
  }
  parlc <- 1-((1-p)^k)
  parlc
}
class(Lc.pareto) <- "theorLc"

Lc.singh <- function(p, parameter=c(2,2))
{
  if(!(is.na(parameter[2]))&(parameter[1]>0)&(parameter[1]<2+parameter[2]))
  {
    b <- parameter[1]-1
    d <- 1/(parameter[2]+1)
  }
  else
  {
    warning("inadmissible parameter. default parameter=c(2,2) is used.")
    b <- 1
    d <- 1/3
  }
  smlc <- pbeta((1-(1-p)^b), (1+d), (b-d))
  smlc
}
class(Lc.singh) <- "theorLc"

Lc.dagum <- function(p, parameter=c(2,2))
{
  if(!(is.na(parameter[2]))&(parameter[1]>1))
  {
    a <- 1/parameter[1]
    b <- parameter[2]
  }
  else
  {
    warning("inadmissible parameter. default parameter=c(2,2) is used.")
    a <- 0.5
    b <- 2
  }
  daglc <- pbeta((p^b), (a+1/b), (1-a))
  daglc
}
class(Lc.dagum) <- "theorLc"

Lc.mehran <- function(x,n)
{
  Lc.min <- Lc(x,n=n)
  p <- Lc.min$p
  L <- Lc.min$L
  k <- length(p)
  q <- c(0,rep(1,k))
  K <- c(rep(0,k),1)
  for(i in k:2)
  {
    q[i] <- 2*p[i]-q[i+1]
  }
  for(i in 2:k)
  {
    K[i] <- 2*L[i-1] - K[i-1]
  }
  beta1 <- (L[2:k]-L[1:(k-1)])/(p[2:k]-p[1:(k-1)])
  beta2 <- (K[2:k]-K[1:(k-1)])/(q[2:k]-q[1:(k-1)])
  beta2 <- beta2[2:(k-1)]
  beta <- rep(0,(k-2))
  for(i in 1:(k-2))
  {
    if(beta1[i]>beta2[i])
      beta[i] <- beta1[i]
    else
    if(beta2[i]>beta1[i+1])
      beta[i] <- beta1[i+1]
    else
      beta[i] <- beta2[i]
  }
  d <- L[2:(k-1)] - beta*p[2:(k-1)]
  if(k==3)
    q <- NULL
  else
    q <- (d[2:(k-2)]-d[1:(k-3)])/(beta[1:(k-3)]-beta[2:(k-2)])
  q <- c(q,1)
  K <- beta*q + d
  L <- c(0,0,K,1)
  p <- c(0,-d[1]/beta[1],q,1)
  L <- L[is.finite(p)]
  p <- p[is.finite(p)]
  Lc.max <- list(p,L)
  names(Lc.max) <- c("p", "L")
  class(Lc.max) <- "Lc"
  Lc.max
}

Lasym <- function(x, n = rep(1, length(x)), interval = FALSE, na.rm = TRUE)
{
    if(!na.rm && any(is.na(x))) return(rep.int(NA_real_, 1L + as.integer(interval)))
    x <- as.numeric(na.omit(x))
    o <- order(x)
    x <- x[o]
    w <- n[o]

    mu <- weighted.mean(x, w)
    xlow <- x < mu
    m <- sum(w[xlow])
    n <- sum(w)
    Lm <- sum(w[xlow] * x[xlow])
    Ln <- sum(w * x)

    if(any(xeq <- x == mu)) {
      a <- sum(w[xeq])
      Lma <- sum(w[xlow | xeq] * x[xlow | xeq])
      Lac <- c(m/n + Lm/Ln, (m + a)/n + Lma/Ln)
      if(!interval) Lac <- mean(Lac)
    } else {
      xm <- max(x[xlow])
      xm1 <- min(x[!xlow])
      delta <- (mu - xm) / (xm1 - xm)
      Lac <- (m + delta)/n + (Lm + delta * xm1)/Ln
    }
    
    Lac
}
