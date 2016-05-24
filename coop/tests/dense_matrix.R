library(coop)
m <- 10
n <- 3

x <- matrix(rnorm(m*n), m, n)
tx <- t(x)
cpx <- crossprod(x)
tcpx <- tcrossprod(x)

y <- x
colnames(y) <- letters[1:ncol(x)]
rownames(y) <- sample(letters, size=nrow(x), replace=TRUE)
ty <- t(y)

check <- function(x, fun1, fun2)
{
  t1 <- fun1(x)
  t2 <- fun2(x)
  stopifnot(all.equal(t1, t2))
}



### Covariance
check(x, cov, coop::covar)
check(tx, cov, coop::covar)
check(tcpx, cov, coop::covar)
check(tcpx, cov, coop::covar)

check(y, cov, coop::covar)
check(ty, cov, coop::covar)



### Pearson correlation
check(x, cor, coop::pcor)
check(tx, cor, coop::pcor)
check(cpx, cor, coop::pcor)
check(tcpx, cor, coop::pcor)

check(y, cor, coop::pcor)
check(ty, cor, coop::pcor)



### Cosine similarity
cosine <- function(x)
{
  cp <- crossprod(x)
  rtdg <- sqrt(diag(cp))
  cos <- cp / tcrossprod(rtdg)
  return(cos)
}

check(x, cosine, coop::cosine)
check(tx, cosine, coop::cosine)
check(cpx, cosine, coop::cosine)
check(tcpx, cosine, coop::cosine)

check(y, cosine, coop::cosine)
check(ty, cosine, coop::cosine)
