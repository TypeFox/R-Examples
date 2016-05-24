## logic for removing/modifying bounds:
## (1) unbounded opt. will have limits of -Inf/Inf
##   [or missing()]
## (2) bounded opt
## fix length mismatch errors!
k <- 19
N <- 20

uniboundtest <- function() {
  m1 <- mle2(k~dbinom(size=N,prob=p),
             start=list(p=0.5))
  m1b <- mle2(k~dbinom(size=N,prob=p),
             start=list(p=0.5),method="L-BFGS-B",upper=0.999)
  p1 <- profile(m1)
  p1b <- profile(m1b)
}
