

# ---------------------------------------------------------------------------------- #

Lmoments <- function(x) {

  camp <- sort(x)
  n <- length(camp)

  nn <- rep(n-1,n)
  pp <- seq(0,n-1)
  p1 <- pp/nn
  p2 <- p1 * (pp-1)/(nn-1)
  p3 <- p2 * (pp-2)/(nn-2)

  b0 <- sum(camp)/n
  b1 <- sum(p1*camp)/n
  b2 <- sum(p2*camp)/n
  b3 <- sum(p3*camp)/n

  l1 <- b0
  l2 <- 2*b1-b0
  lcv <- 2*b1/b0-1
  lca <- 2*(3*b2-b0)/(2*b1-b0)-3
  lkur <- 5*(2*(2*b3-3*b2)+b0)/(2*b1-b0)+6

  Lmom <- c(l1,l2,lcv,lca,lkur)
  names(Lmom) <- c("l1","l2","lcv","lca","lkur")

  return(Lmom)
}


# ---------------------------------------------------------------------------------- #

regionalLmoments <- function(x,cod) {

  # x = vettore contenente i dati di tutte le stazioni
  # cod = vettore  associato ad x con gli identificativi delle stazioni

  if (length(x)!=length(cod)) {stop('x and cod must have the same length')}

  fac <- factor(cod)
  ni <- tapply(x,fac,length)
  k <- nlevels(fac)

  ll <- sapply(split(x,fac),Lmoments)

  l1R <- sum(ni*ll[1,])/sum(ni)
  l2R <- sum(ni*ll[2,])/sum(ni)
  tR <- sum(ni*ll[3,])/sum(ni)
  t3R <- sum(ni*ll[4,])/sum(ni)
  t4R <- sum(ni*ll[5,])/sum(ni)

  rLmom <- c(l1R,l2R,tR,t3R,t4R)
  names(rLmom) <- c("l1R","l2R","lcvR","lcaR","lkurR")

  return(rLmom)
}


# -------------------------------------------------------------------------------- #

LCV <- function (x) {

  # INPUT
  # x = vettore

  x <- sort(x)
  n <- length(x)

  nn <- rep(n-1,n)
  pp <- seq(0,n-1)
  p1 <- pp/nn
  #p2 <- p1 * (pp-1)/(nn-1)
  #p3 <- p2 * (pp-2)/(nn-2)

  b0 <- sum(x)/n
  b1 <- sum(p1*x)/n
  #b2 <- sum(p2*x)/n
  #b3 <- sum(p3*x)/n

  lcv <- 2*b1/b0-1
  #lca <- 2*(3*b2-b0)/(2*b1-b0)-3
  #lkur <- 5*(2*(2*b3-3*b2)+b0)/(2*b1-b0)+6

  return(lcv)
}


# -------------------------------------------------------------------------------- #

LCA <- function (x) {

  # INPUT
  # x = vettore

  x <- sort(x)
  n <- length(x)

  nn <- rep(n-1,n)
  pp <- seq(0,n-1)
  p1 <- pp/nn
  p2 <- p1 * (pp-1)/(nn-1)
  #p3 <- p2 * (pp-2)/(nn-2)

  b0 <- sum(x)/n
  b1 <- sum(p1*x)/n
  b2 <- sum(p2*x)/n
  #b3 <- sum(p3*x)/n

  #lcv <- 2*b1/b0-1
  lca <- 2*(3*b2-b0)/(2*b1-b0)-3
  #lkur <- 5*(2*(2*b3-3*b2)+b0)/(2*b1-b0)+6

  return(lca)
}


# -------------------------------------------------------------------------------- #

Lkur <- function (x) {

  # INPUT
  # x = vettore

  x <- sort(x)
  n <- length(x)

  nn <- rep(n-1,n)
  pp <- seq(0,n-1)
  p1 <- pp/nn
  p2 <- p1 * (pp-1)/(nn-1)
  p3 <- p2 * (pp-2)/(nn-2)

  b0 <- sum(x)/n
  b1 <- sum(p1*x)/n
  b2 <- sum(p2*x)/n
  b3 <- sum(p3*x)/n

  #lcv <- 2*b1/b0-1
  #lca <- 2*(3*b2-b0)/(2*b1-b0)-3
  lkur <- 5*(2*(2*b3-3*b2)+b0)/(2*b1-b0)+6

  return(lkur)
}


