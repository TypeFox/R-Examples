## From file CovSest.R in rrcov

rho.rk2 <- function(x, p, alpha) {
##  x <- t
  z  <- qchisq(1-alpha, p)
  g <- min(z/p - 1, 1)
  uu1 <- (x <= 1-g)
  uu2 <- (x > 1-g & x <= 1+g)
  uu3 <- (x > 1+g)
  zz <- x
  x1 <- x[uu2]
  dd <- ((x1-1)/(4*g))*(3-((x1-1)/g)^2)+.5
  zz[uu1] <- 0
  zz[uu2] <- dd
  zz[uu3] <- 1
  return(zz)
}

rho.rk2.f <- function(x, p, alpha) {
  dq  <- qchisq(1-alpha, p)
  y <- .Fortran("srockech",
    as.double(x),
    as.integer(length(x)),
    as.integer(p),
    as.double(dq),
    PACKAGE="robustvarComp")[[1]]
  return(y)
}

rho.rk2.mean <- function(x, p, alpha) {
  dq  <- qchisq(1-alpha, p)
  y <- .Fortran("drockech",
    as.double(x),
    as.integer(length(x)),
    as.integer(p),
    as.double(dq),
    y=double(1),            
    PACKAGE="robustvarComp")$y
  return(y)
}

fespro.f <- function(sig, fp, fdis, falpha, fdelta) {
  z <- fdis/sig
  dq  <- qchisq(1-falpha, fp)    
  y <- .Fortran("drockech",
    as.double(z),
    as.integer(length(z)),
    as.integer(fp),
    as.double(dq),
    double(1),
    PACKAGE="robustvarComp")[[5]]
    sfun <- y - fdelta
    return(sfun)
}

w.rk2 <- function(x, p, alpha) {
  z <- qchisq(1-alpha, p)
  g <- min(z/p - 1, 1)
  uu1 <- (x <= (1-g) )
  uu2 <- (x > 1-g & x <= 1+g)
  uu3 <- (x > (1+g))
  zz <- x
  x1 <- x[uu2]
  dd <- (3/(4*g))*(1-((x1-1)/g)^2)
  zz[uu1] <- 0
  zz[uu2] <- dd
  zz[uu3] <- 0
  return(zz)
}
