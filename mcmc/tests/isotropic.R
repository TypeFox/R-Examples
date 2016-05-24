library(mcmc)
isotropic <- mcmc:::isotropic
isotropic.logjacobian <- mcmc:::isotropic.logjacobian

# create identity test function
identity <- function(x) x
d.identity <- function(x) 1

# check that isotropic is length preserving for vectors of lengths 1--1000
all(sapply(1:1000, function(x) length(isotropic(identity)(rep(1, x))) == x))
    
# test that isotropic(identity) is an identity function
all.equal(isotropic(identity)(1:10), 1:10)
x <- seq(0, 1, length.out=200)
all.equal(isotropic(identity)(x), x)

# make sure that isotropic.logjacobian(identity, d.identity) is a 0 function
all.equal(isotropic.logjacobian(identity, d.identity)(1:10), 0)

# make sure that 0 as an input does not cause divide-by-zero errors
all.equal(isotropic(identity)(0), 0)
all.equal(isotropic(identity)(0 * 1:4), rep(0, 4))
all.equal(isotropic.logjacobian(identity, d.identity)(0), 0)
all.equal(isotropic.logjacobian(identity, d.identity)(0 * 1:4), 0)

# try isotropic with f(x) = x^2, then we should get 
# istropic(f)(x) := |x| * x
f <- function(x) x^2
all.equal(isotropic(f)(1), 1)
all.equal(isotropic(f)(c(1, 1)), sqrt(2) * c(1, 1))
all.equal(isotropic(f)(c(1, 0, 1)), sqrt(2) * c(1, 0, 1))

# make sure lazy-loading works properly.
g <- function(x) x^2
g.iso <- isotropic(g)
g <- function(x) x
all.equal(g.iso(2), 2*2)
