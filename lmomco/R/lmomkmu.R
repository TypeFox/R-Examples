"lmomkmu" <-
function(para, nmom=5, paracheck=TRUE, tol=1E-6, maxn=100) {

    if(paracheck && ! are.parkmu.valid(para)) {
       warning("Parameters are invalid")
       return()
    }

    names(para$para) <- NULL
    K <- para$para[1]
    M <- para$para[2]

    "marcumq" <- function(a, b, nu=1) {
       pchisq(b^2, df=2*nu, ncp=a^2, lower.tail=FALSE) }

    A <- sqrt(2*K*M)
    B <- sqrt(2*(1+K)*M)
    afunc <- function(x, r=0) {
       Q <- sapply(1:length(x), function(i) { return(marcumq(A,B*x[i], nu=M)) } )
       pdf <- pdfkmu(x, para=para, paracheck=FALSE)
       pdf[is.na(pdf)] <- 0 # This a reverse hack.  The pdf,cdf,qua functions are
       # aggressive on returning NAs for values near or outside the apparent
       # domain of the support, let us consider those zero.
       W <- Q^r * x * pdf
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
    z$source <- "lmomkmu"
    return(z)
}


