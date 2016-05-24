"paremu" <-
function(lmom, checklmom=TRUE, checkbounds=TRUE,
         alsofitT3=FALSE, alsofitT3T4=FALSE, alsofitT3T4T5=FALSE,
         justfitT3T4=TRUE, boundary.tolerance=0.001,
         verbose=FALSE, trackoptim=TRUE) {

   if(alsofitT3T4T5) alsofitT3T4 <- alsofitT3 <- FALSE
   if(alsofitT3T4) alsofitT3   <- FALSE
   if(justfitT3T4) alsofitT3 <- alsofitT3T4 <- alsofitT3T45 <- FALSE

   LARGE <- sqrt(.Machine$double.xmax)


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

   L1 <- lmom$lambdas[1]; L2 <- lmom$lambdas[2]; T2 <- L2/L1
   T3 <- lmom$ratios[3];  T4 <- lmom$ratios[4]; T5 <- lmom$ratios[5]
   lmr <- c(L1, L2, T3, T4)
   if(any(is.na(lmr))) {
      warning("One of the first four L-moments is missing")
      return()
   }
   para <- c(NA, NA)
   names(para) <- c("eta", "mu")
   z <- list(type="emu", para=para, ifail=0, ifailbounds=0,
             message="successful fit", boundarymessage="",
             absboundaryerror = 0, optim=NA,
             source="paremu")

   "mfunc" <- function(e,m) { return(2*m/(1+e^2)) }

   # upper boundary condition
   "t4UP" <- function(t3) {
               return(0.1229 - 0.03548*t3 - 0.1835*t3^2 +
                      2.524*t3^3 - 2.954*t3^4 + 2.001*t3^5 - 0.4746*t3^6)}

   # lower boundary condition
   "t4DN" <- function(t3) {
                return(0.1227  - 0.004433*t3 - 2.845*t3^2 + 18.41*t3^3 - 50.08*t3^4 +
                       83.14*t3^5 - 81.38*t3^6 + 43.24*t3^7 - 9.600*t3^8) }

   if(verbose) message("Testing L-skew and L-kurtosis boundary conditions: ")
   t4u <- t4UP(T3)
   t4d <- t4DN(T3)

   if(checkbounds) {
   if(T4 > t4u) {
      if(verbose) message("   testing whether L-kurtosis, ",T4," is above the upper boundary, ",t4u)
      z$absboundaryerror <- abs(T4 - t4u)
      if(z$absboundaryerror > boundary.tolerance) {
         warning("{T3,T4} lay above the upper boundary upper Eta-Mu")
         z$ifailbounds <- 1
      } else {
         warning("{T3,T4} lay very close to the upper boundary of the Eta-Mu")
         z$ifailbounds <- 0
      }
   }

   if(T4 < t4d) {
      if(verbose) message("   testing whether L-kurtosis, ",T4," is below the lower boundary, ",t4d)
      z$absboundaryerror <- abs(T4 - t4d)
      if(z$absboundaryerror > boundary.tolerance) {
         warning("{T3,T4} lay below the lower boundary upper Eta-Mu")
         z$ifailbounds <- 1
      } else {
         warning("{T3,T4} lay very close to the lower boundary of the Eta-Mu")
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
       e <- exp(par[1]); mu <- exp(par[2]); m <- mfunc(e,mu)
       lmr <- lmomemu(list(para = c(e,mu), type="emu"), nmom=4)
       if(length(lmr) == 0 | is.nan(lmr$ratios[3])) {
          if(trackoptim) message("justfitT3T4: L-moment/parameter blowup")
          return(LARGE)
       }
       SS <- (lmr$ratios[3] - T3)^2 + (lmr$ratios[4] - T4)^2
       SS <- transSS(SS)
       if(trackoptim) {
           message("justfitT3T4: eta = ", round(e,  digits=4),
                   "  mu = ",             round(mu, digits=4),
                   "  m = ",              round(m,  digits=4),
                   "  SS = ",             round(SS, digits=9))
       }
       return(SS)
   }
   "alsofitT3T4T5func" <- function(par) {
       e <- exp(par[1]); mu <- exp(par[2]); m <- mfunc(e,mu)
       lmr <- lmomemu(list(para = c(e,mu), type="emu"), nmom=5)
       if(length(lmr) == 0 | is.nan(lmr$ratios[3])) {
          if(trackoptim) message("alsofitT3T4T5: L-moment/parameter blowup")
          return(LARGE)
       }
       SS <- (lmr$lambdas[1] - L1)^2 + (lmr$lambdas[2]/lmr$lambdas[1] - T2)^2 +
             (lmr$ratios[3]  - T3)^2 + (lmr$ratios[4]  - T4)^2 + (lmr$ratios[5]  - T5)^2
       SS <- transSS(SS)
       if(trackoptim) {
           message("alsofitT3T4T5: eta = ", round(e,  digits=4),
                   "  mu = ",               round(mu, digits=4),
                   "  m = ",                round(m,  digits=4),
                   "  SS = ",               round(SS, digits=9))
       }
       return(SS)
   }
   "alsofitT3T4func" <- function(par) {
       e <- exp(par[1]); mu <- exp(par[2]); m <- mfunc(e,mu)
       lmr <- lmomemu(list(para = c(e,mu), type="emu"), nmom=4)
       if(length(lmr) == 0 | is.nan(lmr$ratios[3])) {
          if(trackoptim) message("alsofitT3T4: L-moment/parameter blowup")
          return(LARGE)
       }
       SS <- (lmr$lambdas[1] - L1)^2 + (lmr$lambdas[2]/lmr$lambdas[1] - T2)^2 +
             (lmr$ratios[3]  - T3)^2 + (lmr$ratios[4]  - T4)^2
       SS <- transSS(SS)
       if(trackoptim) {
           message("alsofitT3T4: eta = ", round(e,  digits=4),
                   "  mu = ",             round(mu, digits=4),
                   "  m = ",              round(m,  digits=4),
                   "  SS = ",             round(SS, digits=9))
       }
       return(SS)
   }
   "alsofitT3func" <- function(par) {
       e <- exp(par[1]); mu <- exp(par[2]); m <- mfunc(e,mu)
       lmr <- lmomemu(list(para = c(e,mu), type="emu"), nmom=3)
       if(length(lmr) == 0 | is.nan(lmr$ratios[3])) {
          if(trackoptim) message("alsofitT3: L-moment/parameter blowup")
          return(LARGE)
       }
       SS <- (lmr$lambdas[1] - L1)^2 + (lmr$lambdas[2]/lmr$lambdas[1] - T2)^2 +
             (lmr$ratios[3]  - T3)^2
       SS <- transSS(SS)
       if(trackoptim) {
           message("alsofitT3: eta = ",   round(e,  digits=4),
                   "  mu = ",             round(mu, digits=4),
                   "  m = ",              round(m,  digits=4),
                   "  SS = ",             round(SS, digits=9))
       }
       return(SS)
   }
   "justfitL1T2func" <- function(par) {
       e <- exp(par[1]); mu <- exp(par[2]); m <- mfunc(e,mu)
       lmr <- lmomemu(list(para = c(e,mu), type="emu"), nmom=3)
       if(length(lmr) == 0 | is.nan(lmr$ratios[3])) {
          if(trackoptim) message("justfitL1T2: L-moment/parameter blowup")
          return(LARGE)
       }
       SS <- (lmr$lambdas[1] - L1)^2 + (lmr$lambdas[2]/lmr$lambdas[1] - T2)^2
       SS <- transSS(SS)
       if(trackoptim) {
           message("justfitL1T2: eta = ", round(e,  digits=4),
                   "  mu = ",             round(mu, digits=4),
                   "  m = ",              round(m,  digits=4),
                   "  SS = ",             round(SS, digits=9))
       }
       return(SS)
   }


   # In the next sequence of steps, we extract the L-moment and parameter table for the
   # Eta-Mu distribution and determine the range of the two parameters. These ranges
   # will be used as lower and upper boundaries in the optimization scheme. Next,
   # the Euclidean error term is determined according to three logical arguments to the
   # function. Finally, the initial value or guesses are extracted.
   EMU <- .lmomcohash$EMU_lmompara_byeta[complete.cases(.lmomcohash$EMU_lmompara_byeta),]
   ETA.min <- min(EMU$ETA);   ETA.max <- max(EMU$ETA)
   MU.min  <- min(EMU$MU);    MU.max  <- max(EMU$MU)
   if(verbose) message("Tabulated eta limits: ",ETA.min," and ",ETA.max)
   if(verbose) message("Tabulated  mu limits: ",MU.min," and ",MU.max)
   err <- 0
   if(justfitT3T4) {
       err <- (EMU$T3 - T3)^2 + (EMU$T4 - T4)^2
       afunc <- justfitT3T4func
   } else if(alsofitT3T4T5) {
       err <- (EMU$L1 - L1)^2 + (EMU$L2/EMU$L1 - T2)^2 + (EMU$T3 - T3)^2 + (EMU$T4 - T4)^2 + (EMU$T5 - T5)^2
       afunc <- alsofitT3T4T5func
   } else if(alsofitT3T4) {
       err <- (EMU$L1 - L1)^2 + (EMU$L2/EMU$L1 - T2)^2 + (EMU$T3 - T3)^2 + (EMU$T4 - T4)^2
       afunc <- alsofitT3T4func
   } else if(alsofitT3) {
       err <- (EMU$L1 - L1)^2 + (EMU$L2/EMU$L1 - T2)^2 + (EMU$T3 - T3)^2
       afunc <- alsofitT3func
   } else {
       err <- (EMU$L1 - L1)^2 + (EMU$L2/EMU$L1 - T2)^2
       afunc <- justfitL1T2func
   }
   #err <- (EMU$L1 - L1)^2 + (EMU$L2/EMU$L1 - T2)^2 + (EMU$T3 - T3)^2
   ix  <- 1:length(EMU[,1])
   ixe <- ix[err == min(err, na.rm=TRUE)]
   ETA.guess <- EMU$ETA[ixe]
   MU.guess  <- EMU$MU[ixe]
   if(verbose) message("Initial eta = ", ETA.guess, " and mu = ",MU.guess)

   rt <- NULL
   try( rt <- optim(c(log(ETA.guess), log(MU.guess)), afunc))#,
                    #method="L-BFGS-B",
                    #lower=log(c(ETA.min,MU.min)),
                    #upper=log(c(ETA.max,MU.max))) )
   if(is.null(rt)) {
      z$ifail <- 10
      z$message <- "optimizer returned NULL, failure"
      return(z)
   } else {
      z$optim <- rt
   }

   ETA <- exp(rt$par[1])
   MU  <- exp(rt$par[2])

   para <- c(ETA,MU)
   z$para <- para

   M <- mfunc(ETA,MU)
   if(MU < M/2 | MU > M) {
      z$ifail <- 12
      warning("Parameter M is outside of the proper range for the Eta-Mu distribution")
      return(z)
   }

   return(z)
}

