dinvgauss <- function(x, mu = stop("no shape arg"), lambda = 1)
{
#  Density of inverse Gaussian distribution
#  GKS  15 Jan 98
#
        if(any(mu<=0)) stop("mu must be positive")
        if(any(lambda<=0)) stop("lambda must be positive")
        d <- ifelse(x>0,sqrt(lambda/(2*pi*x^3))*exp(-lambda*(x-mu)^2/(2*mu^2*x)),0)
        if(!is.null(Names <- names(x)))
                names(d) <- rep(Names, length = length(d))
        d
}

pinvgauss <- function(q, mu = stop("no shape arg"), lambda = 1)
{
#  Inverse Gaussian distribution function
#  GKS  15 Jan 98
#
        if(any(mu<=0)) stop("mu must be positive")
        if(any(lambda<=0)) stop("lambda must be positive")
        n <- length(q)
        if(length(mu)>1 && length(mu)!=n) mu <- rep(mu,length=n)
        if(length(lambda)>1 && length(lambda)!=n) lambda <- rep(lambda,length=n)
        lq <- sqrt(lambda/q)
        qm <- q/mu
        p <- ifelse(q>0,pnorm(lq*(qm-1))+exp(2*lambda/mu)*pnorm(-lq*(qm+1)),0)
        if(!is.null(Names <- names(q)))
                names(p) <- rep(Names, length = length(p))
        p
}

rinvgauss <- function(n, mu = stop("no shape arg"), lambda = 1)
{
#  Random variates from inverse Gaussian distribution
#  Reference:
#      Chhikara and Folks, The Inverse Gaussian Distribution,
#      Marcel Dekker, 1989, page 53.
#  GKS  15 Jan 98
#
        if(any(mu<=0)) stop("mu must be positive")
        if(any(lambda<=0)) stop("lambda must be positive")
        if(length(n)>1) n <- length(n)
        if(length(mu)>1 && length(mu)!=n) mu <- rep(mu,length=n)
        if(length(lambda)>1 && length(lambda)!=n) lambda <- rep(lambda,length=n)
        y2 <- rchisq(n,1)
        u <- runif(n)
        r1 <- mu/(2*lambda) * (2*lambda + mu*y2 - sqrt(4*lambda*mu*y2 + mu^2*y2^2))
        r2 <- mu^2/r1
        ifelse(u < mu/(mu+r1), r1, r2)
}

qinvgauss  <- function(p, mu = stop("no mean arg"), lambda = 1)
{
#  Quantiles of the inverse Gaussian distribution
#  Dr Paul Bagshaw
#  Centre National d'Etudes des Telecommunications (DIH/DIPS)
#  Technopole Anticipa, France
#  paul.bagshaw@cnet.francetelecom.fr
#  23 Dec 98
#
  if(any(mu <= 0.))
    stop("mu must be positive")
  if(any(lambda <= 0.))
    stop("lambda must be positive")
  n <- length(p)
  if(length(mu) > 1 && length(mu) != n)
    mu <- rep(mu, length = n)
  if(length(lambda) > 1 && length(lambda) != n)
    lambda <- rep(lambda, length = n)
  thi <- lambda / mu
  U <- qnorm (p)
  r1 <- 1 + U / sqrt (thi) + U^2 / (2 * thi) + U^3 / (8 * thi * sqrt(thi))
  x <- r1
  for (i in 1:10) {
    cum <- pinvgauss (x, 1., thi)
    dx <- (cum - p) / dinvgauss (x, 1., thi)
    dx <- ifelse (is.finite(dx), dx, ifelse (p > cum, -1, 1))
    dx[dx < -1] <- -1
    if (all(dx == 0.)) break
    x <- x - dx
  }
  x * mu
}
