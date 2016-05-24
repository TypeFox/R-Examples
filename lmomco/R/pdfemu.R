"pdfemu" <-
function(x, para, paracheck=TRUE) {
   if(paracheck == TRUE) {
      if(! are.paremu.valid(para)) return()
   }

   SMALL <- 1E-6

   E <- para$para[1]
   M <- para$para[2]
   h <- 1 / (1-E^2)
   H <- E / (1-E^2)
   nu <- M - 1/2
   nu.is.negative.integer <- FALSE
   if(nu < 0) {
     if(as.integer(nu) == nu) nu.is.negative.integer <- TRUE
   }

   f <- sapply(1:length(x), function(i) {
          xi  <- x[i]
          if(is.na(xi) || xi < 0) return(NA)
          if(! is.finite(xi)) return(0)
          xx <- xi^2
          if(E <= SMALL) { # eta going to zero
             return(2*(2*M)^(2*M) / gamma(2*M) * xi^(2*(2*M) - 1) * exp(-2*M*xx))
          } else {
             #message("Other ")
             toI <- 2*M*H*xx
             tmpB <- 4*sqrt(pi)*M^(M+1/2)*h^M  / ( gamma(M)*H^nu )
             tmpB <- tmpB * xi^(2*M) * exp(-2*M*h*xx)
             B <- besselI(toI, nu=nu)
             #message("X      ",xi)
             #message("toI    ",toI)
             #message("tmpB   ",tmpB)
             #message("B      ",B,"\n")
             if(is.finite(B) & tmpB != 0) { # eta going to +/- unity
                return(tmpB * B)
             } else {
                return(2*M^M / gamma(M) * xi^(2*M-1) * exp(-M*xx))
             }
          } })
   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}

