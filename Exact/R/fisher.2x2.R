fisher.2x2<-function(x, alternative="two.sided"){
  m <- sum(x[, 1])
  n <- sum(x[, 2])
  k <- sum(x[1, ])
  x <- x[1, 1]
  PVAL <- switch(alternative,less = phyper(x, m, n, k),
                 greater = phyper(x - 1, m, n, k, lower.tail = FALSE),
                 two.sided = {
                   relErr <- 1 + 10^(-7)
                   lo <- max(0, k - n)
                   hi <- min(k, m)
                   support <- lo:hi
                   d <- dhyper(support, m, n, k, log = TRUE)
                   d <- exp(d - max(d))
                   d <- d/sum(d)
                   sum(d[d <= d[x - lo + 1] * relErr])
                 })
  PVAL
}
