cosine <- function(x, y)
{
  cp <- crossprod(x, y)
  cp / sqrt(crossprod(x) * crossprod(y))
}


x <- rnorm(30)
y <- rnorm(30)

t1 <- as.vector(cosine(x, y))
t2 <- coop::cosine(x, y)
stopifnot(all.equal(t1, t2, check.attributes=FALSE))


t1 <- cor(x, y)
t2 <- coop::pcor(x, y)
stopifnot(all.equal(t1, t2, check.attributes=FALSE))


t1 <- cov(x, y)
t2 <- coop::covar(x, y)
stopifnot(all.equal(t1, t2, check.attributes=FALSE))
