"lmomst3" <-
function(para, bypoly=TRUE) {
    if(! are.parst3.valid(para)) return()

    U <- para$para[1]
    A <- para$para[2]
    N <- para$para[3]

    SMALL.NU <- 1.000001 # arrived from manual experiments
    LARGE.NU <- 1000     # limits of experiments yielding the polynomial
    if(N < SMALL.NU) N <- SMALL.NU
    if(N > LARGE.NU) N <- LARGE.NU

    L1   <- U
    L2   <- 2^(6-4*N)*pi*A*sqrt(N)*exp(lgamma(2*N - 2) - 4*lgamma(N/2))
    TAU3 <- 0
    TAU4 <- NA
    TAU5 <- 0
    "poly6" <- function(t4) {
    	 b <- 0.001599337
    	ce <- c(  0.085903249,  2.749217116, -6.723564788, 13.445446378,
    	        -17.717402013, 14.586006328, -6.788165377,  1.360977224)
    	t6 <- b + ce[1]*t4^1 + ce[2]*t4^2 + ce[3]*t4^3 + ce[4]*t4^4 +
                  ce[5]*t4^5 + ce[6]*t4^6 + ce[7]*t4^7 + ce[8]*t4^8
        return(t6)
    }
    
    if(bypoly) {
       "polyt4" <- function(nu) {
          lnu <- log(nu)
            b <- -6.962256e-04
          ce <- c(-1.812284e+00,  6.299075e-01, -6.859132e-02, -2.266309e-02,
                   9.921283e-03, -1.651713e-03,  1.342971e-04, -4.399008e-06)
          lgt4 <- b + ce[1]*lnu^1 + ce[2]*lnu^2 + ce[3]*lnu^3 + ce[4]*lnu^4 +
                      ce[5]*lnu^5 + ce[6]*lnu^6 + ce[7]*lnu^7 + ce[8]*lnu^8
          tau4 <- exp(lgt4)
          return(tau4)
       }
       TAU4 <- polyt4(N)
    } else if(N >= LARGE.NU) { # treat as a normal distribution
       TAU4 <- 30/pi * atan(sqrt(2)) - 9
    } else {
       TAU4 <- (15/2)*exp(lgamma(N) - lgamma(1/2) - lgamma(N - 1/2))
       afunc <- function(x) {
          AA <- 1/sqrt(x) * (1 - x)^(N - (3/2)) * pbeta(x, 0.5, N/2)^2
          return(AA)
       }
       int <- NULL
       try(int <- integrate(afunc, 0, 1), silent=TRUE)
       if(is.null(int)) {
          warning("Bad integration")
       } else {
          TAU4 <- TAU4 * int$value - (3/2)
          if(TAU4 >= 1) TAU4 <- NA
       }
    }
    TAU6 <- poly6(TAU4)
    lam <- c(L1, L2,    TAU3*L2, TAU4*L2, TAU5*L2, TAU6*L2)
    rat <- c(NA, L2/L1, TAU3, TAU4, TAU5, TAU6)
    names(rat) <- NULL
    names(lam) <- NULL
    zz <- list(lambdas=lam, ratios=rat,
               trim=0, leftrim=NULL, rightrim=NULL,
               source="lmomst3")
    return(zz)
}

