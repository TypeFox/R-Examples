mh<-function (tab, Gamma = 1)
{
  tab1<-tab
  d1 <- dim(tab1)
  stopifnot((d1[1] == 2) & (d1[2] == 2))
  if (length(d1) == 2)
    tab1 <- array(tab1, c(2, 2, 1))
  gamma <- Gamma
  stopifnot(sum(as.vector(tab1))>0)
  tab1<-tab1[,,apply(tab1,3,sum)>0]
  if (length(dim(tab1))==2) tab1<-array(tab1,c(2,2,1))

  check2x2xktable<-function(tab){
    r1<-tab[1,1,]+tab[1,2,]
    r2<-tab[2,1,]+tab[2,2,]
    c1<-tab[1,1,]+tab[2,1,]
    c2<-tab[1,2,]+tab[2,2,]
    mc<-pmin(pmin(r1,r2),pmin(c1,c2))
    if (max(mc)==0) {
      warning("One of the 2x2 or 2x2xk tables is degenerate.")
      stop("There is a problem with your data.")
    }
  }
  check2x2xktable(tab1)


  one2x2 <- function(tb, gamma = gamma) {
    m1 <- tb[1, 1] + tb[1, 2]
    m2 <- tb[2, 1] + tb[2, 2]
    n <- tb[1, 1] + tb[2, 1]
    mx <- min(n, m1)
    mn <- max(0, n - m2)
    x <- mn:mx
    g <- rep(0, mx + 1)
    pr <- dFNCHypergeo(x, m1, m2, n, gamma)
    g[(mn + 1):(mx + 1)] <- pr
    g
  }
  one2x2xk <- function(tbk, gamma = gamma) {
    k <- dim(tbk)[3]
    g <- 1
    for (i in 1:k) {
      gi <- one2x2(tbk[, , i], gamma = gamma)
      g <- gconv(g, gi)
    }
    g
  }
  g1 = one2x2xk(tab1, gamma = gamma)
  val1 <- (0:(length(g1) - 1))
  actuala <- sum(tab1[1, 1, ])
  pa <- sum(g1[val1 >= actuala])
  list(pval=pa,A=actuala)
}
