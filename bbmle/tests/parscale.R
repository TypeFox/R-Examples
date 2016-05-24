library(bbmle)
old_opt <- options(digits=3)
tracelevel <- 0

## source("~/lib/R/pkgs/bbmle/pkg/R/mle.R

set.seed(1002)
X <- rexp(1000, rate = 0.0001)
f <- function(X, rate) {
  if (tracelevel>0 && rate<0) cat("rate<0: ",rate,"\n")
  -sum(dexp(X, rate = rate, log = TRUE))
}
if (FALSE) {
  ## L-BFGS-B violates bounds, and gets stuck at lower bound
  m <- mle2(minuslogl = f,
            data = list(X = X),
            start = list(rate = 0.01),
            method = "L-BFGS-B",
            control = list(trace = tracelevel, 
                           parscale = 1e-4),
            lower = c(rate = 1e-9))

  profile(m, std.err=0.0001) ## finds new optimum

  fsc <- function(X, rate) {
    -sum(dexp(X, rate = rate*1e-4, log = TRUE))
  }
  msc <- mle2(minuslogl = fsc,
            data = list(X = X),
            start = list(rate = 100),
            method = "L-BFGS-B",
            control = list(trace = tracelevel),
            lower = c(rate = 1e-5))

  ## does it work if we scale by hand?
  ##   no, identical problem
}

## works fine with a better starting point
m <- mle2(minuslogl = f,
          data = list(X = X),
          start = list(rate = 0.001),
          method = "L-BFGS-B",
          control = list(trace = tracelevel,
                         parscale=1e-4),
              lower = c(rate = 1e-9))
vcov(m)
confint(m)


## works OK despite warnings about 1-dimensional opt. with N-M
(m0 <- mle2(minuslogl = f,
          data = list(X = X),
          start = list(rate = 0.01),
          method = "Nelder-Mead",
          control = list(trace = tracelevel, parscale = 1e-4)))
vcov(m0)

confint(m0)
confint(m0,method="quad")
## very similar (good quadratic surface, not surprising)

m1 <- mle2(minuslogl = f,
          data = list(X = X),
          start = list(rate = 0.01),
          method = "BFGS",
          control = list(trace = tracelevel, parscale = 1e-4))


## gets stuck? will have to investigate ...
m2 <- mle2(minuslogl = f,
           data = list(X = X),
           start = list(rate = 0.01),
           optimizer = "optimize",
           lower=1e-9,upper=0.1)

vcov(m2)
options(old_opt)
