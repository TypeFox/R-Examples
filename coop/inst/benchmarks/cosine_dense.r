cosine_R <- function(x)
{
  cp <- crossprod(x)
  rtdg <- sqrt(diag(cp))
  cos <- cp / tcrossprod(rtdg)
  return(cos)
}

library(compiler)
cosine_R <- cmpfun(cosine_R)

library(coop)
library(rbenchmark)
cols <- cols <- c("test", "replications", "elapsed", "relative")
reps <- 25

m <- 10000
n <- 250
x <- matrix(rnorm(m*n), m, n)

benchmark(cosine_R(x), cosine(x), replications=reps, columns=cols)
### Too slow!
#benchmark(lsa::cosine(x), cosine(x), replications=reps, columns=cols)
