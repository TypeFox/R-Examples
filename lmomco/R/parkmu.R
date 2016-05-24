"parkmu" <-
function(lmom, checklmom=TRUE, checkbounds=TRUE,
         alsofitT3=FALSE, alsofitT3T4=FALSE, alsofitT3T4T5=FALSE,
         justfitT3T4=TRUE, boundary.tolerance=0.001,
         verbose=FALSE, trackoptim=TRUE) {

   if(alsofitT3T4T5) alsofitT3T4 <- alsofitT3 <- FALSE
   if(alsofitT3T4)   alsofitT3   <- FALSE
   if(justfitT3T4)   alsofitT3 <- alsofitT3T4 <- alsofitT3T45 <- FALSE

   if(length(lmom$L1) == 1) {  # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
   }
   if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
   }

   if(length(lmom$lambdas) < 5) {
      warning("Five or more L-moments are required even if optional settings do not use all five")
   }

   LARGE <- sqrt(.Machine$double.xmax)

   L1 <- lmom$lambdas[1]; L2 <- lmom$lambdas[2]; T2 <- L2/L1
   T3 <- lmom$ratios[3];  T4 <- lmom$ratios[4]; T5 <- lmom$ratios[5]
   lmr <- c(L1, L2, T3, T4)
   if(any(is.na(lmr))) {
      warning("One of the first four L-moments is missing")
      return()
   }
   para <- c(NA, NA)
   names(para) <- c("kappa", "mu")
   z <- list(type="kmu", para=para, ifail=0, ifailbounds=0,
             message="successful fit", boundarymessage="",
             absboundaryerror = 0, diracdelta=0, optim=NA,
             source="parkmu")

   "mfunc" <- function(k,m) { return(m*(1+k)^2/(1+2*k)) }

   # upper boundary condition
   "t4UP" <- function(t3) {
                return(0.1227  - 0.004433*t3 - 2.845*t3^2 + 18.41*t3^3 - 50.08*t3^4 +
                       83.14*t3^5 - 81.38*t3^6 + 43.24*t3^7 - 9.600*t3^8) }

   "t4DNA" <- function(t3) {
                 return(0.1226 - 0.3206*t3 - 102.4*t3^2 - 4.753E4*t3^3 +
                        -7.605E6*t3^4 - 5.244E8*t3^5 - 1.336E10*t3^6)}

   "t4DNB" <- function(t3) {
                 return(0.09328  - 1.488*t3 + 16.29*t3^2 - 205.4*t3^3 +
                        1545*t3^4 - 5595*t3^5 + 7726*t3^6)}

   "t4DNC" <- function(t3) {
                 return(0.07245 - 0.8631*t3 + 2.031*t3^2 - 0.01952*t3^3 +
                        -0.7532*t3^4 + 0.7093*t3^5 - 0.2156*t3^6)}

   if(verbose) message("Testing L-skew and L-kurtosis boundary conditions: ")
   t4u  <- t4UP(T3)
   t4dA <- t4DNA(T3)
   t4dB <- t4DNB(T3)
   t4dC <- t4DNC(T3)

   if(checkbounds) {
   if(T4 > t4u) {
      if(verbose) message("   testing whether L-kurtosis, ",T4," is above the upper boundary, ",t4u)
      z$absboundaryerror <- abs(T4 - t4u)
      if(z$absboundaryerror > boundary.tolerance) {
         warning("{T3,T4} lay above the upper boundary upper Kappa-Mu")
         z$ifailbounds <- 1
      } else {
         #warning("{T3,T4} lay very close to the upper boundary of the Kappa-Mu")
         z$ifailbounds <- 0
      }
   }
   if((-0.0169 <= T3 & T3 <= 0) & (0.1218 <= T4 & T4 <= 0.1271) & T4 < t4dA) {
      if(verbose) message("   testing whether L-kurtosis ",T4," is below the lower boundary A, ",t4dA)
      if(abs(T4 - t4dA) > boundary.tolerance) {
         warning("{T3,T4} lay outside the left-most lower boundary of the Kappa-Mu")
         z$ifailbounds <- 2
      } else {
         #warning("{T3,T4} lay very close to the left-most lower boundary of the Kappa-Mu")
         z$ifailbounds <- 0
      }
   }
   if((-0.0169 <= T3 & T3 <= 0.2041) & (-0.0226 <= T4 & T4 <= 0.1218) & T4 < t4dB) {
      if(verbose) message("   testing whether L-kurtosis ",T4," is below the lower boundary B, ",t4dB)
      if(abs(T4 - t4dB) > boundary.tolerance) {
         z$ifailbounds <- 3
         warning("{T3,T4} lay outside the lower-left lower boundary of the Kappa-Mu")
      } else {
         #warning("{T3,T4} lay very close to the lower-left lower boundary of the Kappa-Mu")
         z$ifailbounds <- 0
      }
   }
   if((0.2041 <= T3 & T3 <= 1) & (-0.0226 <= T4 & T4 <= 1) & T4 < t4dC) {
      if(verbose) message("   testing whether L-kurtosis ",T4," is below the lower boundary C, ",t4dC)
       if(abs(T4 - t4dC) > boundary.tolerance) {
          z$ifailbounds <- 4
          warning("{T3,T4} lay outside the lower-right lower boundary of the Kappa-Mu")
       } else {
          #warning("{T3,T4} lay very close to the lower-right lower boundary of the Kappa-Mu")
          z$ifailbounds <- 0
       }
   }

   if(z$ifailbounds > 0) {
      z$ifail <- 1
      z$message <- "fit rejected"
      return(z)
   } else {
      if(verbose) message("   boundary conditions seem satisfied")
   }
   }

   "transSS" <- function(SS) return(sqrt(SS))

   # The following functions to be optimized over
   "justfitT3T4func" <- function(par) {
       k <- exp(par[1]); mu <- exp(par[2]); m <- mfunc(k,mu)
       lmr <- lmomkmu(list(para = c(k,mu), type="kmu"), nmom=4)
       if(length(lmr) == 0 | is.nan(lmr$ratios[3])) {
          if(trackoptim) message("justfitT3T4: L-moment/parameter blowup")
          return(LARGE)
       }
       SS <- (lmr$ratios[3] - T3)^2 + (lmr$ratios[4] - T4)^2
       SS <- transSS(SS)
       if(trackoptim) {
           message("justfitT3T4: kappa = ", round(k,  digits=4),
                   "  mu = ",               round(mu, digits=4),
                   "  m = ",                round(m,  digits=4),
                   "  SS = ",               round(SS, digits=9))
       }
       return(SS)
   }
   "alsofitT3T4T5func" <- function(par) {
       k <- exp(par[1]); mu <- exp(par[2]); m <- mfunc(k,mu)
       lmr <- lmomkmu(list(para = c(k,mu), type="kmu"), nmom=5)
       if(length(lmr) == 0 | is.nan(lmr$ratios[3])) {
          if(trackoptim) message("alsofitT3T4T5: L-moment/parameter blowup")
          return(LARGE)
       }
       SS <- (lmr$lambdas[1] - L1)^2 + (lmr$lambdas[2]/lmr$lambdas[1] - T2)^2 +
             (lmr$ratios[3]  - T3)^2 + (lmr$ratios[4]  - T4)^2 + (lmr$ratios[5]  - T5)^2
       SS <- transSS(SS)
       if(trackoptim) {
           message("alsofitT3T4T5: kappa = ", round(k,  digits=4),
                   "  mu = ",                 round(mu, digits=4),
                   "  m = ",                  round(m,  digits=4),
                   "  SS = ",                 round(SS, digits=9))
       }
       return(SS)
   }
   "alsofitT3T4func" <- function(par) {
       k <- exp(par[1]); mu <- exp(par[2]); m <- mfunc(k,mu)
       lmr <- lmomkmu(list(para = c(k,mu), type="kmu"), nmom=4)
       if(length(lmr) == 0 | is.nan(lmr$ratios[3])) {
          if(trackoptim) message("alsofitT3T4: L-moment/parameter blowup")
          return(LARGE)
       }
       SS <- (lmr$lambdas[1] - L1)^2 + (lmr$lambdas[2]/lmr$lambdas[1] - T2)^2 +
             (lmr$ratios[3]  - T3)^2 + (lmr$ratios[4]  - T4)^2
       SS <- transSS(SS)
       if(trackoptim) {
           message("alsofitT3T4: kappa = ", round(k,  digits=4),
                   "  mu = ",               round(mu, digits=4),
                   "  m = ",                round(m,  digits=4),
                   "  SS = ",               round(SS, digits=9))
       }
       return(SS)
   }
   "alsofitT3func" <- function(par) {
       k <- exp(par[1]); mu <- exp(par[2]); m <- mfunc(k,mu)
       lmr <- lmomkmu(list(para = c(k,mu), type="kmu"), nmom=3)
       if(length(lmr) == 0 | is.nan(lmr$ratios[3])) {
          if(trackoptim) message("alsofitT3: L-moment/parameter blowup")
          return(LARGE)
       }
       SS <- (lmr$lambdas[1] - L1)^2 + (lmr$lambdas[2]/lmr$lambdas[1] - T2)^2 +
             (lmr$ratios[3]  - T3)^2
       SS <- transSS(SS)
       if(trackoptim) {
           message("alsofitT3: kappa = ", round(k,  digits=4),
                   "  mu = ",             round(mu, digits=4),
                   "  m = ",              round(m,  digits=4),
                   "  SS = ",             round(SS, digits=9))
       }
       return(SS)
   }
   "justfitL1T2func" <- function(par) {
       k <- exp(par[1]); mu <- exp(par[2]); m <- mfunc(k,mu)
       lmr <- lmomkmu(list(para = c(k,mu), type="kmu"), nmom=3)
       if(length(lmr) == 0 | is.nan(lmr$ratios[3])) {
          if(trackoptim) message("justfitL1T2: L-moment/parameter blowup")
          return(LARGE)
       }
       SS <- (lmr$lambdas[1] - L1)^2 + (lmr$lambdas[2]/lmr$lambdas[1] - T2)^2
       SS <- transSS(SS)
       if(trackoptim) {
           message("justfitL1T2: kappa = ", round(k,  digits=4),
                   "  mu = ",               round(mu, digits=4),
                   "  m = ",                round(m,  digits=4),
                   "  SS = ",               round(SS, digits=9))
       }
       return(SS)
   }


   # In the next sequence of steps, we extract the L-moment and parameter table for the
   # Kappa-Mu distribution and determine the range of the two parameters. These ranges
   # will be used as lower and upper boundaries in the optimization scheme. Next,
   # the Euclidean error term is determined according to three logical arguments to the
   # function. Finally, the initial value or guesses are extracted.
   KMU <- .lmomcohash$KMU_lmompara_bykappa[complete.cases(.lmomcohash$KMU_lmompara_bykappa),]
   KAPPA.min <- min(KMU$KAPPA);    KAPPA.max <- max(KMU$KAPPA)
   MU.min    <- min(KMU$MU);       MU.max    <- max(KMU$MU)
   if(verbose) message("Tabulated kappa limits: ",KAPPA.min," and ",KAPPA.max)
   if(verbose) message("Tabulated    mu limits: ",MU.min," and ",MU.max)
   err <- 0
   if(justfitT3T4) {
       err <- (KMU$T3 - T3)^2 + (KMU$T4 - T4)^2
       afunc <- justfitT3T4func
   } else if(alsofitT3T4T5) {
       err <- (KMU$L1 - L1)^2 + (KMU$L2/KMU$L1 - T2)^2 + (KMU$T3 - T3)^2 + (KMU$T4 - T4)^2 + (KMU$T5 - T5)^2
       afunc <- alsofitT3T4T5func
   } else if(alsofitT3T4) {
       err <- (KMU$L1 - L1)^2 + (KMU$L2/KMU$L1 - T2)^2 + (KMU$T3 - T3)^2 + (KMU$T4 - T4)^2
       afunc <- alsofitT3T4func
   } else if(alsofitT3) {
       err <- (KMU$L1 - L1)^2 + (KMU$L2/KMU$L1 - T2)^2 + (KMU$T3 - T3)^2
       afunc <- alsofitT3func
   } else {
       err <- (KMU$L1 - L1)^2 + (KMU$L2/KMU$L1 - T2)^2
       afunc <- justfitL1T2func
   }
   #err <- (KMU$L1 - L1)^2 + (KMU$L2/KMU$L1 - T2)^2 + (KMU$T3 - T3)^2
   ix  <- 1:length(KMU[,1])
   ixe <- ix[err == min(err, na.rm=TRUE)]
   KAPPA.guess <- KMU$KAPPA[ixe]
   MU.guess    <- KMU$MU[ixe]
   if(verbose) message("Initial kappa = ", KAPPA.guess, " and mu = ",MU.guess)

   rt <- NULL
   try( rt <- optim(log(c(KAPPA.guess, MU.guess)), afunc))#,
                    #method="L-BFGS-B",
                    #lower=log(c(KAPPA.min,MU.min)),
                    #upper=log(c(KAPPA.max,MU.max))))
   if(is.null(rt)) {
      z$ifail <- 10
      z$message <- "optimizer returned NULL, failure"
      return(z)
   } else {
      z$optim <- rt
   }

   KAPPA <- exp(rt$par[1])
   MU    <- exp(rt$par[2])

   para <- c(KAPPA,MU)
   z$para <- para
   z$diracdelta <- pdfkmu(0, z)

   M <- MU*(1+KAPPA)^2/(1+2*KAPPA)
   if(MU < 0 | MU > M) {
      z$ifail <- 12
      warning("Parameter M is outside of the proper range for the Kappa-Mu distribution")
      return(z)
   }

   return(z)
}


