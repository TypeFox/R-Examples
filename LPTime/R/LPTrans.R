
LPTrans <- function(x,m = 4){
    u <- ecdf(x)(x)
    leg.f <- legendre.polynomials(m, normalized=FALSE)
    X <- polynomial.values(leg.f, (2*u - 1))
    S <- scale(do.call(cbind,X)[,-1])
}

