"pargep" <-
function(lmom, checklmom=TRUE, checkdomain=TRUE, maxit=10, verbose=FALSE) {

   "canGEPsingle" <- function(t2, t3) {
      "T3geplo" <- function(t2) {
          b  <- 0.17395362
          me <- c(-0.07186818, 0.13659663, 2.91130745, -4.81008588, 3.49997828, -0.8016708)
          t3 <- b + me[1]*t2^1 + me[2]*t2^2 + me[3]*t2^3 + me[4]*t2^4 + me[5]*t2^5 + me[6]*t2^6
          return(t3)
      }
      "T3gephiA" <- function(t2) {
          b  <- 0.3268947
          me <- c(0.4306951, -3.3853219, 11.5967566, -17.6576974, 14.3702814, -4.7305213)
          t3 <- b + me[1]*t2^1 + me[2]*t2^2 + me[3]*t2^3 + me[4]*t2^4 + me[5]*t2^5 + me[6]*t2^6
          return(t3)
      }
      "T3gephiB" <- function(t2) {
          b  <- 6.825674e+00
          me <- c(-3.034089e+02, 5.548407e+03, -4.976073e+04, 2.195468e+05, -3.820133e+05)
          t3 <- b + me[1]*t2^1 + me[2]*t2^2 + me[3]*t2^3 + me[4]*t2^4 + me[5]*t2^5
          return(t3)
      }
       lo <-  T3geplo(t2)
      hiA <- T3gephiA(t2)
      hiB <- T3gephiB(t2)

      if(t2 > 0.8636418) return(FALSE) # this point is a hacked intersection of the polynomial boundaries
      if(t2 > 0.825)     return(FALSE) # yes another test, but this appears the real edge of numerical ability
      if(t2 < 0.085)     return(FALSE)
      if(t3 <= lo)       return(FALSE)
      if(t2 >= 0.110083) {
         if(t3 >= hiA)   return(FALSE)
      } else {
         if(t3 >= hiB)   return(FALSE)
      }
      return(TRUE)
   }

   "canGEP" <- function(t2, t3) {
      lmr <- list(t2=t2, t3=t3)
       zz <- sapply(1:length(t2), function(i) { return(canGEPsingle(lmr$t2[i], lmr$t3[i]) ) })
      return(zz)
   }


   para <- vector(mode = "numeric", length = 3)
   names(para) <- c("beta", "kappa", "h")

   if(length(lmom$L1) == 1) { # convert to named L-moments
     lmom <- lmorph(lmom)     # nondestructive conversion!
   }
   if(checklmom & ! are.lmom.valid(lmom)) {
     warning("L-moments are invalid")
     return()
   }

   L1 <- lmom$lambdas[1]; L2 <- lmom$lambdas[2]
   L3 <- lmom$lambdas[3]; L4 <- lmom$lambdas[4]
   T2 <-  lmom$ratios[2]; T3 <-  lmom$ratios[3]

   if(L1 <= 0) {
      warning("The mean is not > 0")
      return()
   }

   if(checkdomain) {
      can.do.GEP <- canGEP(T2, T3)
      if(! can.do.GEP) {
         if(verbose) message("The L-moments are outside the numerical domain of the GEP")
         para <- rep(NA, 3)
         return(list(type="gep", para=para, convergence=10, error=NA, its=NA, source="pargep"))
      }
   }

   hatE11 <- L1;             hatE22 <- L2 + hatE11
   hatE33 <- (L3 + 3*hatE22 - hatE11)/2
   #hatE44 <- (L4 + 10*hatE33 - 6*hatE22 + hatE11)/5 # not used
   hatE22oE11 <- hatE22/hatE11
   hatE33oE11 <- hatE33/hatE11
   #hatE44oE11 <- hatE44/hatE11  # not used

   EPS <- 1E-4
   beta.guess <- NULL # experimental feature in case of user knowledge, disabled
   beta.guess <- ifelse(is.null(beta.guess), hatE11, beta.guess)

   # Function to be used for Beta estimation
   "bfunc" <- function(b, k=NULL, h=NULL) {
      gep <- vec2par(c(10^b,k,h), type="gep")
      E11 <- expect.max.ostat(1, para=gep, qua=quagep)
      return(abs(E11 - hatE11))
   }

   "khfunc" <- function(kh, b=NULL) {
      k <- 10^kh[1]; h <- 10^kh[2]
      gep <- vec2par(c(b,k,h), type="gep")
      E33 <- expect.max.ostat(3, para=gep, qua=quagep)
      E22 <- expect.max.ostat(2, para=gep, qua=quagep)
      E11 <- expect.max.ostat(1, para=gep, qua=quagep)
      return(sqrt((E33/E11 - hatE33oE11)^2 + (E22/E11 - hatE22oE11)^2))
   }

   "bkhfunc" <- function(bkh) {
      b <- 10^bkh[1]; k <- 10^bkh[2]; h <- 10^bkh[3]
      gep <- vec2par(c(b,k,h), type="gep")
      E33 <- expect.max.ostat(3, para=gep, qua=quagep)
      E22 <- expect.max.ostat(2, para=gep, qua=quagep)
      E11 <- expect.max.ostat(1, para=gep, qua=quagep)
      return(sqrt((E33 - hatE33)^2 + (E22 - hatE22)^2 + (E11 - hatE11)^2))
   }

   if(! verbose) message("STATUS (iteration): ", appendLF=FALSE)
   khpar <- c(1,1) # Very frequently a solution can be hit with log10{k,h} = {1,1}
   it <- 0 
   while(1) {
      it <- it + 1
      if(! verbose) message(it,"-", appendLF=FALSE)

      tmp <- NULL 
      try(tmp <- optim(khpar, khfunc, control=list(trace=0), b=beta.guess), silent=TRUE)
      khpar <- runif(2, min=-1, max=1) # if we hit the next iteration, have some different start points
      if(is.null(tmp)) tmp <- list(convergence = 9999)
      if(tmp$convergence != 0) { # failed, so let us gradual open up a random starting points so function is not deterministic
         if(verbose) message("KH joint optimization failed, new starting point (B)")
         try(tmp <- optim(runif(2, min=-0.5, max=1.5), khfunc, control=list(trace=0), b=beta.guess), silent=TRUE)
         if(is.null(tmp)) tmp <- list(convergence = 9999)
         if(tmp$convergence != 0) {
            if(verbose) message("KH joint optimization failed, new starting point (C)")
            try(tmp <- optim(runif(2, min=-1, max=2), khfunc, control=list(trace=0), b=beta.guess), silent=TRUE)
            if(is.null(tmp)) tmp <- list(convergence = 9999)
            if(tmp$convergence != 0) {
               if(verbose) message("KH joint optimization failed, new starting point (D)")
               try(tmp <- optim(runif(2, min=-1.5, max=2.5), khfunc, control=list(trace=0), b=beta.guess), silent=TRUE)
               if(is.null(tmp)) tmp <- list(convergence = 9999)
               if(tmp$convergence != 0) {
                  if(verbose) message("KH joint optimization failed, exiting")
                  para <- rep(NA, 3)
                  zz <- list(type="gep", para=para, convergence=2, error=NA, its=it, source="pargep")
                  return(zz)
               }
            }
         }
      }
      #if(verbose) {
      #  print("KH OPTIMIZATION")
      #  print(tmp)
      #}
      K <- 10^tmp$par[1];  H <- 10^tmp$par[2]

      # Now solve for Beta
      tmp <- NULL
      try(tmp <- optimize(bfunc, c(-2,6), k=K, h=H), silent=TRUE)
      if(is.null(tmp)) {
         if(verbose) message("Optimization on E11 using parameter b failed")
         para <- rep(NA, 3)
         zz <- list(type="gep", para=para, convergence=3, error=NA, its=it, source="pargep")
         return(zz)
      }
      #if(verbose) {
      #   print("B OPTIMIZATION")
      #   print(tmp)
      #}
      B <- 10^tmp$minimum
      beta.guess <- B
      para <- c(B, K, H)

      lmrgep <- lmomgep(vec2par(para, type='gep'))
      error <- abs(lmrgep$lambdas[2] - L2)/L2 + abs(lmrgep$lambdas[3] - L3)/L3
      converge <- ifelse(error < EPS, 0, 1)
      if(converge == 0 | it >= maxit) break
   }
   if(converge != 0) {
      if(! verbose) message("adding special last attempt and incrementing iterations")
      it <- it + 1
      tmp <- NULL
      try(tmp <- optim(log10(para), bkhfunc, control=list(trace=0)), silent=TRUE)
      if(is.null(tmp)) tmp <- list(convergence = 9999)
      if(tmp$convergence != 0) {
         if(verbose) message("Optimization on E33, E22, and E11 failed")
         para <- rep(NA, 3)
         zz <- list(type="gep", para=para, convergence=4, error=NA, its=it, source="pargep")
         return(zz)
      }
      B <- 10^tmp$par[1]; K <- 10^tmp$par[2];  H <- 10^tmp$par[3]
      para <- c(B, K, H)
      lmrgep <- lmomgep(vec2par(para, type='gep'))
      error <- abs(lmrgep$lambdas[2] - L2)/L2 + abs(lmrgep$lambdas[3] - L3)/L3
      converge <- ifelse(error < EPS, 0, 1)
   } else {
      if(! verbose) message("done")
   }
   zz <- list(type="gep", para=para, convergence=converge, error=error, its=it, source="pargep")
   return(zz)
}
