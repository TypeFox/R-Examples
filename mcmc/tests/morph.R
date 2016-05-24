library(mcmc)
isotropic <- mcmc:::isotropic
isotropic.logjacobian <- mcmc:::isotropic.logjacobian

# make sure morph identity works properly
TestMorphIdentity <- function(m.id) {
  ident.func <- function(x) x
  if (!all.equal(m.id$transform(1:10), 1:10))
    return(FALSE)
  if (!all.equal(m.id$inverse(1:10), 1:10))
    return(FALSE)
  x <- seq(-1,1, length.out=15)
  if (!all.equal(sapply(x, m.id$lud(function(x) dnorm(x, log=TRUE))),
                 dnorm(x, log=TRUE)))
    return(FALSE)
  if (!all.equal(m.id$outfun(ident.func)(x), x))
    return(FALSE)
  return(TRUE)
}

TestMorphIdentity(morph())
TestMorphIdentity(morph.identity())

TestMorphIdentityOutfun <- function(m) {
  f <- m$outfun(NULL)
  x <- 1:20
  if (!identical(x, f(x)))
    return(FALSE)
  f <- m$outfun(c(6, 8))
  if (!identical(x[c(6, 8)], f(x)))
    return(FALSE)
  i <- rep(FALSE, 20)
  i[c(1, 3, 5)] <- TRUE
  f <- m$outfun(i)
  if (!identical(x[i], f(x)))
    return(FALSE)
  return(TRUE)
}

TestMorphIdentityOutfun(morph())
TestMorphIdentityOutfun(morph.identity())

# make sure that morph and morph.identity give back the same things
all.equal(sort(names(morph.identity())), sort(names(morph(b=1))))

# test center parameter, univariate version
zero.func <- function(x) 0
center <- 2
x <- seq(-1,1, length.out=15)
morph.center <- morph(center=center)
all.equal(sapply(x, morph.center$transform), x-center)
all.equal(sapply(x, morph.center$inverse), x+center)
all.equal(sapply(x, morph.center$lud(function(y) dnorm(y, log=TRUE))),
          dnorm(x, log=TRUE, mean=-2))

# test center parameter, multivariate version
center <- 1:4
x <- rep(0, 4)
morph.center <- morph(center=center)
lud.mult.dnorm <- function(x) prod(dnorm(x, log=TRUE))
all.equal(morph.center$transform(x), x-center)
all.equal(morph.center$inverse(x), x+center)
all.equal(morph.center$lud(lud.mult.dnorm)(x),
          lud.mult.dnorm(x - center))
# test 'r'.
r <- 1
morph.r <- morph(r=r)
x <- seq(-1, 1, length.out=20)
all.equal(sapply(x, morph.r$lud(function(x) dnorm(x, log=TRUE))),
          dnorm(x, log=TRUE))
x <- seq(1.1, 2, length.out=10)
all(sapply(x, morph.r$lud(function(x) dnorm(x, log=TRUE)))
    !=
    dnorm(x, log=TRUE))

TestExponentialEvenPWithRInverse <- function() {
  r <- 0.3
  p <- 2.2
  morph.r <- morph(r=r, p=p)
  x <- seq(0, r, length.out=20)
  all.equal(x, sapply(x, morph.r$inverse))
}

TestExponentialEvenPWithRInverse()

# make sure morph$lud passes '...' arguments.
mean <- 2
ident.morph <- morph()
dnorm.morph <- ident.morph$lud(function(x, mean=0)
                                 dnorm(x, mean=mean, log=TRUE))
all.equal(dnorm.morph(2, mean), dnorm(2, mean=mean, log=TRUE))
x <- seq(-3, 3, length.out=20)
m2 <- morph(r=10)
dnorm.morph <- m2$lud(function(x, mean)
                        dnorm(x, mean=mean, log=TRUE))
all.equal(sapply(x, function(y) dnorm.morph(y, 2)),
          dnorm(x, mean=2, log=TRUE))

# make sure morph$outfun passes '...' arguments.
outfun.orig <- function(x, mean) x + mean
ident.morph <- morph()
mean <- 1
outfun.morph <- ident.morph$outfun(outfun.orig)
all.equal(outfun.morph(1:10, mean), 1:10+mean)

m2 <- morph(r=10)
outfun.morph <- m2$outfun(outfun.orig)
all.equal(sapply(1:10, function(x) outfun.morph(x, mean)), 1:10+mean)

###########################################################################
# test built-in exponential and polynomial transformations.
f <- morph(b=3)
x <- seq(0, 10, length.out=100)
all.equal(x, sapply(sapply(x, f$transform), f$inverse))

f <- morph(p=3)
all.equal(x, sapply(sapply(x, f$transform), f$inverse))

f <- morph(p=3, r=10)
all.equal(-10:10, Vectorize(f$transform)(-10:10))

f <- morph(p=3, b=1)
all.equal(x, sapply(sapply(x, f$transform), f$inverse))
