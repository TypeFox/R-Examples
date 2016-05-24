"pdfkmu" <-
function(x, para, paracheck=TRUE) {
   if(paracheck == TRUE & ! are.parkmu.valid(para)) return()

   SMALL <- .Machine$double.eps

   # x typically is a normalized signal level (r/rrms) where
   # r is the received signal and rrms is the root mean square
   # of the signal (the standard deviation)
   K  <- para$para[1]
   MU <- para$para[2]
   fixedM <- FALSE
   if(! is.finite(K)) { M <- MU; fixedM <- TRUE }
   if(fixedM) {
      tmpA <- 4*M/exp(2*M)
      tmpB <- sqrt(2*M*pi)/exp(M)
      toI <- 4*M*0
      B1 <- besselI(toI, nu=1)
      if(! is.finite(B1)) B1 <- 0
      B2 <- besselI(M, nu=1/2)
      if(! is.finite(B2)) B2 <- 0
      diracdelta <- tmpA*B1 + (1 - tmpB*B2)
      names(diracdelta) <- "Dirac Delta x=0"
      # This resetting of the Dirac Delta functions part of the parameter list is made
      # so that the conditional tests numerical equivalence is effectively bypassed.
      # This feature is provided so that should a user create their own parameter list manually
      # and not through the function vec2par that it will be assumed that the Dirac computed here is fine.
      if(length(para$diracdelta) == 0 | is.na(para$diracdelta)) para$diracdelta <- diracdelta
      if(diracdelta != para$diracdelta) {
          warning("Dirac delta (x=0) computed herein does not match that embedded in the ",
                  "parameter object, going to use the freshly computed one")
          warning("Dirac delta = ", diracdelta,"   and embedded = ", para$diracdelta)
      }
      #message("Note: The Dirac Delta function for (x=0) for this parameterized Kappa-Mu ",
      #"distribution provides ",round(diracdelta, digits=6)," of total probability.\n")
   }

   f <- sapply(1:length(x), function(i) {
                 xi  <- x[i]
                 if(xi < 0)          return(NA)
                 if(! is.finite(xi)) return(NA)
                 if(is.na(xi))       return(NA)
                 if(fixedM) {
                    tmpA <- 4*M/exp(2*M*(1+xi^2))
                    tmpB <- sqrt(2*M*pi)/exp(M)
                    toI <- 4*M*xi
                    ops <- options(warn=-1)
                    B1 <- besselI(toI, nu=1)
                    #if(! is.finite(B1)) B1 <- 0
                    B2 <- besselI(M, nu=1/2)
                    #if(! is.finite(B2)) B2 <- 0
                    options(ops)
                    V <- tmpA*B1 + (1 - tmpB*B2)
                    ifelse(xi < SMALL, return(    diracdelta),
                                       return(V - diracdelta))
                 } else {
                    if(MU == 1) {
                       tmpB <- 2 * (1+K) / exp(K) * xi * exp( -(1+K) * xi^2)
                       toI <- 2 * sqrt(K*(1+K)) * xi
                       ops <- options(warn=-1)
                       B <- besselI(toI, nu=0, expon.scaled=TRUE)/exp(-toI)
                       options(ops)
                    } else if(K == 0) {
                       tmpB <- 2*MU^MU / gamma(MU) * xi^(2*MU - 1) * exp( -MU * xi^2 )
                       B <- 1
                    } else {
                       tmpA <- 2 * MU * ( 1 + K )^( (MU+1) / 2 ) / ( K^( (MU-1) / 2 ) * exp( MU * K ) )
                       tmpB <- tmpA * xi^MU * exp( -MU * (1+K) * xi^2 )
                       toI  <- 2 * MU * sqrt( K * (K+1) ) * xi
                       ops <- options(warn=-1)
                       B <- besselI(toI, nu=MU-1, expon.scaled=TRUE)/exp(-toI)
                       options(ops)
                    }
                 }
                 f <- tmpB * B
                 ifelse(is.finite(f), return(f), return(NA))
              })

   names(f) <- NULL
   f[! is.finite(f)] <- NA
   f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}

