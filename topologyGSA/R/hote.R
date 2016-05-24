.hote <- function(y1, y2, exact, cliques=NULL) {
  y1.num   <- nrow(y1)
  y2.num   <- nrow(y2)
  gene.num <- ncol(y1)

  y1.bar <- colMeans(y1)
  y2.bar <- colMeans(y2)
  y1.s   <- cov(y1)
  y2.s   <- cov(y2)

  y.diff <- y1.bar - y2.bar
  if (!is.null(cliques)) {
    y1.s <- qpIPF(y1.s, cliques)
    y2.s <- qpIPF(y2.s, cliques)
  }

  k <- y1.num + y2.num - gene.num - 1

  if (exact) {
    s  <- ((y1.num-1)*y1.s + (y2.num-1)*y2.s) / (y1.num + y2.num - 2)
    t2 <- ((y1.num*y2.num) / (y1.num+y2.num)) * (y.diff %*% solve(s) %*% y.diff)

    t.obs <- as.vector( t2 * k / (gene.num * (y1.num + y2.num - 2)) )
    alpha.obs <- 1 - pf(t.obs, gene.num, k)

    list(alpha.obs=alpha.obs,
         t.obs=t.obs,
         df=c(gene.num, k))

  } else {

    s <- y1.s/y1.num + y2.s/y2.num
    t2 <- as.vector( ((y1.num*y2.num) / (y1.num+y2.num)) * (y.diff %*% solve(s) %*% y.diff) )

    list(t.obs=t2,
         df=c(gene.num, k))
  }
}

.hotePaired <- function(y1, y2, cli.moral, perm=FALSE) {
  y1.num <- nrow(y1)
  y.diff <- y1 - y2

  if (perm) {
    signs <- matrix(sample(c(1,-1), y1.num*ncol(y1), replace=TRUE),
                    nrow=y1.num)
    y.diff <- y.diff * signs
  }

  y.bar <- colMeans(y.diff)
  y.centr <- y.diff - y.bar
  y.s <- qpIPF(cov(y.diff), cli.moral)
  t2 <- as.numeric(y1.num * (t(y.bar) %*% solve(y.s) %*% y.bar))

  p <- ncol(y1)
  np <- y1.num - p

  t2 * np / (p * (y1.num-1))
}
