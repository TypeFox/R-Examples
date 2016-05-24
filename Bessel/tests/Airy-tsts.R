library(Bessel)

### Replicate some testing "ideas" from  zqcai.f (TOMS 644 test program)

##  Generates airy functions and their derivatives from zairy
##  and zbiry and checks them against the wronskian evaluation in the
##  region -pi/3 <= arg(z) <= pi/3:
##
##              Ai(z)*Bi'(z)-Ai'(z)*Bi(z) = 1/pi.
##
##  in the remainder of the cut plane, the identities
##
##    Ai(z)  = sqrt(-z)*( J(-1/3,zr) + J(1/3,zr) )/3
##
##    Ai'(z) =        z*( J(-2/3,zr) - J(2/3,zr) )/3
##
##    Bi(z)  = i* sqrt(-z/3) *( c1* H(1/3,1,zr) - c2* H(1/3,2,zr) )/2
##
##    Bi'(z) = i*(-z)/sqrt(3)*( c2* H(2/3,1,zr) - c1* H(2/3,2,zr) )/2
##
##  are checked where  zr = (2/3)(-z)^(3/2)  with
##  c1 = exp(pi*i/6),
##  c2 = conjg(c1)    and i**2 = -1.

allEQ <- function(a,b) all.equal(a,b, tol = 1e-12)

### Wronskian  Ai(z)*Bi'(z) - Ai'(z)*Bi(z) == 1/pi  ----------------
N <- 100
I.pi <- rep.int(1/pi + 0*1i, N)
c1 <- exp(pi * 1i/6) ## = sqrt(3)/2 + i/2
c2 <- Conj(c1)

set.seed(101)
for(n in 1:50) {
    cat(".")
    z <- complex(re = rnorm(N),
                 im = rnorm(N))
    ai  <- AiryA(z)
    dai <- AiryA(z, deriv=1)
    bi  <- AiryB(z)
    dbi <- AiryB(z, deriv=1)
    stopifnot(allEQ(ai * dbi - dai * bi, I.pi))

    ## The remaining checks are only valid in this z-plane "sector":
    Lz <- abs(Arg(z)) > pi/3
    z <- z[Lz]
    ai  <-  ai[Lz]
    dai <- dai[Lz]
    bi  <-  bi[Lz]
    dbi <- dbi[Lz]

    zr <- 2/3 * (-z)^(3/2)

##    Ai(z)  = sqrt(-z)*( J(-1/3,zr) + J(1/3,zr) )/3
##    Ai'(z) =        z*( J(-2/3,zr) - J(2/3,zr) )/3

    if(FALSE) ## FIXME -- needs BesselJ() for *NEGATIVE* nu :
    stopifnot(allEQ(ai, sqrt(-z)*(BesselJ(zr, -1/3) + BesselJ(zr, 1/3))/3),
              allEQ(dai,       z*(BesselJ(zr, -2/3) - BesselJ(zr, 2/3))/3))

##    Bi(z)  = i* sqrt(-z/3) *( c1* H(1/3,1,zr) - c2* H(1/3,2,zr) )/2
##    Bi'(z) = i*(-z)/sqrt(3)*( c2* H(2/3,1,zr) - c1* H(2/3,2,zr) )/2

    stopifnot(allEQ(bi, 1i*sqrt(-z/3)*
                    (c1*BesselH(1, zr, 1/3) - c2*BesselH(2, zr, 1/3))/2),
              allEQ(dbi, -1i*z/sqrt(3)*
                    (c2*BesselH(1, zr, 2/3) - c1*BesselH(2, zr, 2/3))/2))

}; cat("\n")
