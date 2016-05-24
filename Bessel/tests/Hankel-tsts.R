library(Bessel)

### Replicate some testing "ideas" from  zqcbh.f (TOMS 644 test program)

##  zqcbh generates sequences of H Bessel functions for kind 2 from
##  zbesh and checks them against analytic continuation formulas
##  in the (z,fnu) space:
##
##  kode = 1 tests (analytic continuation formulae, i**2 = -1):
##
##  H(fnu,2,z)= -exp(i*pi*fnu) * H(fnu,1,-z),    -pi < arg(z) <= 0
##
##            = 2*cos(pi*fnu) * H(fnu,2,-z) + exp(i*pi*fnu)*H(fnu,1,-z),
##
##                                                0 < arg(z) <= pi
##  kode = 2 tests for kinds 1 and 2:
##
##         exp(-i*z) * H(fnu,1,z) = [exp(-i*z)* H(fnu,1,z)]
##
##         exp( i*z) * H(fnu,2,z) = [exp( i*z)* H(fnu,2,z)]
##
##  where the left side is computed with kode = 1 and the right side
##  with kode = 2.

N <- 100 # speed ...
set.seed(1)
for(nu in round(sort(rlnorm(20)), 2)) {
    cat("nu =", format(nu),"  ")
    for(n in 1:5) {
        cat(".")
        z <- complex(re = rnorm(N),
                     im = rnorm(N))
        b1  <- BesselH(1, z, nu=nu)
        b1. <- BesselH(1, z, nu=nu, expon.scaled = TRUE)
        b2  <- BesselH(2, z, nu=nu)
        b2. <- BesselH(2, z, nu=nu, expon.scaled = TRUE)
        exp.i.pn <- exp(1i*pi*nu)
        b2.alt <-
            ifelse(0 < Arg(z),
                   2*cos(pi*nu) * BesselH(2,-z,nu) + exp.i.pn * BesselH(1,-z,nu),
                   ## else Arg(z) <= 0:
                   -exp.i.pn * BesselH(1,-z,nu))
        stopifnot(all.equal(b2,  b2.alt,            tol = 30e-16),
                  all.equal(exp(-1i * z) * b1, b1., tol = 80e-16),
                  all.equal(exp( 1i * z) * b2, b2., tol = 80e-16))
    }
    cat("\n")
}


cat('Time elapsed: ', proc.time(),'\n') # "stats"

if(!interactive()) warnings()
