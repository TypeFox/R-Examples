"cdfemu" <-
function(x, para, paracheck=TRUE, yacoubsintegral=TRUE) {
    if(paracheck == TRUE) {
      if(! are.paremu.valid(para)) return()
    }

    ETA.SMALL <- 1E-3
    SMALL <- 1E-10
    ONEminusSMALL <- 1 - SMALL
    LARGE <- sqrt(.Machine$double.xmax)

    names(para$para) <- NULL
    E <- para$para[1]
    M <- para$para[2]
    b <- sqrt( 2*M / (1-E^2) )
    Ymuleft <- sqrt(pi)*2^(3/2 - M)*(1 - E^2)^M
    Ymuleft <- Ymuleft / (gamma(M) * E^(M-1/2))

    if(E >= 1 - ETA.SMALL) {
       #warning("Eta is near unity, going to use integration of the pdf")
       yacoubsintegral <- FALSE
    }
    if(E <= ETA.SMALL) {
       #warning("Eta is near zero, going to use integration of the pdf")
       yacoubsintegral <- FALSE
    }
    "Ymu" <- function(x=NULL, y=NULL) {
       if(y < 0) return(NA)
       "afunc" <- function(t) {
           B <- vector(mode="numeric", length=length(t))
           for(i in 1:length(B)) {
              toI <- t[i]^2*x
              b <- NULL
              try(b <- besselI(toI, nu=M-1/2, expon.scaled=TRUE)/exp(-toI),
                               silent=TRUE)
              if(is.null(b) | is.nan(b) | ! is.finite(b)) b <- LARGE
              B[i] <- b
           }
           B[! is.finite(B)] <- LARGE
           z <- exp(-t^2)*t^(2*M) * B
           z[! is.finite(z)] <- LARGE
           return( z )
       }
       int1 <- NULL
       try( int1 <- integrate(afunc, lower=y, upper=Inf), silent=TRUE )
       if(is.null(int1)) return(NA)
       val <- int1$value
       if(val == 0) return(0)
       return(Ymuleft * int1$value)
    }

    f <- sapply(1:length(x), function(i) {
             xi <- x[i]
             if(xi == 0) return(0)
             if(is.na(xi) || xi < 0) return(NA) # return zero?
             if(yacoubsintegral) {
                Y <- Ymu(x=E, y=b*xi)
                if(is.na(Y)) {# | Y <= SMALL) {
                   return(ONEminusSMALL)
                } else {
                   return(1 - Y)
                }
             } else {
                int1 <- NULL
                try( int1 <- integrate(pdfemu, 0, xi, para=para, paracheck=FALSE) )
                if(is.null(int1)) return(NA)
                val <- int1$value
                ifelse(val >= ONEminusSMALL, return(ONEminusSMALL),
                                             return(val))  } })
    names(f) <- NULL
    return(f)
}

