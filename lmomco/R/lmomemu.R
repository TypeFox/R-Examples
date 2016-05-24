"lmomemu" <-
function(para, nmom=5, paracheck=TRUE, tol=1E-6, maxn=100) {

    if(paracheck && ! are.paremu.valid(para)) {
       warning("Parameters are invalid")
       return()
    }

    names(para$para) <- NULL
    E <- para$para[1]
    M <- para$para[2]


    b <- sqrt( 2*M / (1-E^2) )
    Ymuleft <- sqrt(pi)*2^(3/2 - M)*(1 - E^2)^M
    Ymuleft <- Ymuleft / (gamma(M) * E^(M-1/2))
    "yacoubsintegral" <- function(x=NULL, y=NULL) {
       if(y < 0) return(NA)
       "afunc" <- function(t) {
           B <- vector(mode="numeric", length=length(t))
           for(i in 1:length(B)) {
              toI <- t[i]^2*x
              B[i] <- besselI(toI, nu=M-1/2, expon.scaled=TRUE)/exp(-toI)
           }
           B[! is.finite(B)] <- .Machine$double.xmax
           z <- exp(-t^2)*t^(2*M) * B
           return( z )
       }
       int1 <- NULL
       try( int1 <- integrate(afunc, lower=y, upper=Inf), silent=TRUE )
       if(is.null(int1)) return(NA)
       return( Ymuleft * int1$value )
    }


    h <- 1 / (1 - E^2)
    A <- E
    B <- sqrt(2*h*M)
    afunc <- function(x, r=0) {
       Y <- sapply(1:length(x), function(i) { return(yacoubsintegral(A,B*x[i])) } )
       pdf <- pdfemu(x, para=para, paracheck=FALSE)
       pdf[is.na(pdf)] <- 0 # This a reverse hack.  The pdf,cdf,qua functions are
       # aggressive on returning NAs for values near or outside the apparent
       # domain of the support, let us consider those zero.
       W <- Y^r * x * pdf
       return(W)
    }

    alphas <- rep(NA, length(nmom))
    for(r in 1:nmom) {
         int <- NULL
         try(int <- integrate(afunc, 0, Inf, subdivisions=100, r=(r-1)), silent=TRUE)
         if(is.null(int)) {
            alphas[r] <- NaN;
            next;
         }
         alpha <- int$value
         alphas[r] <- alpha
    }
    z <- pwm2lmom(pwm.alpha2beta(alphas))
    z$source <- "lmomemu"
    return(z)
}

