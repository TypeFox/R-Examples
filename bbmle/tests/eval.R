## I am experiencing difficulties with one of my modeling function (bbmle::mle2)
## which, like other modeling functions in R, uses match.call() to
## retrieve and save the original function call for future use.
## I'll describe the problem for bbmle and then show that I can
## provoke a similar problem with lm().

## ============
## PART I: mle2()

  library(bbmle)

  x <- 0:10
  y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
  d <- data.frame(x,y)

## The key is to call the modeling function from within another
## function which passes additional arguments via ... 

  ff <- function(d,...) {
    mle2(y~dpois(lambda=ymean),start=list(ymean=mean(y)),data=d,...)
  }

  ff(d)
  try(ff(d,control=list(maxit=1000)))

##   Error in call$control$parscale : 
##      object of type 'symbol' is not subsettable

## This happens when I try:

##  call$control$parscale <- eval.parent(call$control$parscale)

## in 'normal' circumstances call$control and call$control$parscale
## are either NULL or well-specified ...

## Debugging mle2 shows that the results of match.call() are

##   mle2(minuslogl = y ~ dpois(lambda = ymean), start = list(ymean = mean(y)), 
##       data = d, control = ..1)

## ============
## PART II: lm()

## I can find a similar issue with lm(), although admittedly
## I have to work a bit harder/do something a little bit more
## obscure.

  L1 <-  lm(y~1,data=d,tol=1e-6)
  L1$call

  ff2 <- function(d,...) {
    lm(y~1,data=d,...)
  }

  tt <- 1e-6
  L2 <- ff2(d,tol=tt)
  L2$call

  try(update(L2,.~.+x))

## Error in lm.fit(x, y, offset = offset, singular.ok = singular.ok, ...) : 
##   ..1 used in an incorrect context, no ... to look in

 ## similar issue in curve3d().  How does curve() work?


