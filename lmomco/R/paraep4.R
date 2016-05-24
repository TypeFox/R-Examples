"paraep4" <-
function(lmom, checklmom=TRUE,
         method=c("A", "DG", "ADG"),
         sqrt.t3t4=TRUE, eps=1e-4,
         checkbounds=TRUE, kapapproved=TRUE,
         A.guess=NULL, K.guess=NULL, H.guess=NULL) {

    method <- match.arg(method)

    if(checklmom && ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid.")
      return()
    }

    if(length(lmom$L1) != 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }

    para <- vector(mode="numeric", length=4)
    names(para) <- c("xi","alpha","kappa","h")
    para <- c(NA, NA, NA, NA);


    z <- list(type   = 'aep4',    para   = para,
              source = "paraep4", method = method,
              ifail  = 0, ifailtext="")

    L234 <- list(para_L234 = para,
                 ifail_L234=NA,
                 optim.converg_L234=NA,
                 optim.message_L234=NA,
                 optim.value_L234=NA,
                 optim.counts_L234=NA)
    T34 <- list(para_T34 = para,
                ifail_T34=NA,
                optim.converg_T34=NA,
                optim.message_T34=NA,
                optim.value_T34=NA,
                optim.counts_T34=NA)

    L1 <- lmom$lambdas[1]
    L2 <- lmom$lambdas[2]
    L3 <- lmom$lambdas[3]
    L4 <- lmom$lambdas[4]
    T2 <- lmom$ratios[2]
    T3 <- lmom$ratios[3]
    T4 <- lmom$ratios[4]


    if(checkbounds) {
       a <- abs(T3)
       co <- c(0.7755464,  -3.3354852,  14.1955782, -29.9090294,
               37.2141451, -24.7411869,   6.7997646)
       T4.lowerbounds <-  co[1]*a   + co[2]*a^2 + co[3]*a^3 +
                          co[4]*a^4 + co[5]*a^5 + co[6]*a^6 + co[7]*a^7

       if(T4 < T4.lowerbounds) {
          if(kapapproved) {
             z <- parkap(lmom)
             z$ifailkap  <- z$ifail
             z$ifailtextkap <- z$ifailtext
             z$ifailtext <- "TAU4 is estimated as below limits of AEP4, Kappa fit instead"
             z$source <- "paraep4 --> parkap"
             return(z)
          } else {
             z$ifail     <- 4
             z$ifailtext <- "TAU4 is estimated as below limits of AEP4, Kappa not fit instead"
             return(z)
          }
       }
    }

    if(is.null(A.guess)) A.guess <- 1

    err <- (T3 - .lmomcohash$AEPkh2lmrTable$T3)^2 +
           (T4 - .lmomcohash$AEPkh2lmrTable$T4)^2
    if(is.null(K.guess)) {
       K.guess <- .lmomcohash$AEPkh2lmrTable$K[err == min(err)]
    }
    if(is.null(H.guess)) {
       H.guess <- .lmomcohash$AEPkh2lmrTable$H[err == min(err)]
    }

    para.guess <- vec2par(c(0,A.guess,K.guess,H.guess), type="aep4")
    if(! are.paraep4.valid(para.guess)) {
       message <- "One or more of the guesses of A, K, and H (regardless of method choice) are invalid."
       warning(message)
       z$ifail <- 3
       z$ifailtext <- message
       return(z)
    }

    if(method != "A") {
      opt <- NULL
        "fn" <- function(ps, ...) {
             para <- list(para=c(0,ps), type="aep4")
             slmr <- lmomaep4(para, paracheck=FALSE)
             return(log(1 + (L2 - slmr$lambdas[2])^2
                          + (L3 - slmr$lambdas[3])^2
                          + (L4 - slmr$lambdas[4])^2))
        }
        try( opt <- optim(c(A.guess,K.guess,H.guess), fn), silent=TRUE)
          if(is.null(opt) | length(opt$par) == 0) {
            L234$ifail_L234 <- 1
            L234$optim.converg_L234 <- NA
            L234$optim.message_L234 <- NA
            L234$optim.value_L234   <- NA
            L234$optim.counts_L234  <- NA
            z$L234 <- L234
            if(method == "DG") {
               message <- "The function optim failed or reports failure on its own behalf."
               warning(message)
               z$ifail <- 1
               z$ifailtext <- message
               return(z)
            }
          } else {
            para[2:4] <- opt$par
            A <- para[2]
            K <- para[3]
            H <- para[4]
            KmK <- 1/K - K
            H2H1 <- exp(lgamma(2/H) - lgamma(1/H))
            U <- L1 - A * KmK * H2H1
            para[1] <- U

            L234$para_L234 <- para
            L234$ifail_L234 <- opt$convergence
            if(method == "DG") {
                z$para <- para
                z$ifail <- L234$ifail_L234
            }
            L234$optim.converg_L234 <- opt$convergence
            L234$optim.message_L234 <- opt$message
            L234$optim.value_L234   <- opt$value
            L234$optim.counts_L234  <- opt$counts
            z$L234 <- L234

           if(method == "DG") {
              if(! are.paraep4.valid(z)) {
                 message <- "One or more parameters are not valid: Delicado-Goria method."
                 warning(message)
                 z$ifailtext <- message
                 z$ifail <- 3
              }
              return(z)
           }
         }
      }

      "sqrtit" <- function(x) { return(x) }
      if(sqrt.t3t4) "sqrtit" <- function(x) { return(sqrt(x)) }

       opt <- NULL
      "fn" <- function(ps, ...) {
           para <- list(para=c(0,1,ps), type="aep4")
           #print(para)
           slmr <- lmomaep4(para, paracheck=FALSE, t3t4only=TRUE)
           return(sqrtit((T3 - slmr$T3)^2 + (T4 - slmr$T4)^2))
      }
      try( opt <- optim(c(K.guess,H.guess), fn), silent=TRUE)
         if(is.null(opt) | length(opt$par) == 0) {
            T34$ifail_T34 <- 1
            T34$optim.converg_T34 <- NA
            T34$optim.message_T34 <- NA
            T34$optim.value_T34   <- NA
            T34$optim.counts_T34  <- NA
            z$T34 <- T34
            message <- "The function optim failed or reports failure on its own behalf."
            warning(message)
            z$ifail <- 1
            z$ifailtext <- message
            return(z)
         }
         para[3:4] <- opt$par
         K <- para[3]
         H <- para[4]
         KmK <- 1/K - K
         H2H1 <- exp(lgamma(2/H) - lgamma(1/H))
         Ihalf <- pbeta(1/2, shape1=1/H, shape2=2/H)
         KK    <- K*K
         KKK   <- KK*K

         L2a <- -K * KmK^2              / (1+KK)
         L2b <-  2 * KK * (1/KKK + KKK) / (1+KK)^2 * Ihalf
         A <- L2 / ((L2a + L2b) * H2H1)
         para[2] <- A
         U <- L1 - A * KmK * H2H1
         para[1] <- U

         z$para <-  para
         T34$para_T34 <- para
         T34$ifail_T34 <- opt$convergence
         z$ifail <- T34$ifail_T34
         if(opt$value > eps) {
            message <- "Judging a solution failure based on eps value to convergence error: one of the A methods."
            warning(message)
            z$ifailtext <- message
            z$ifail <- 2
         }
         T34$optim.converg_T34 <- opt$convergence
         T34$optim.message_T34 <- opt$message
         T34$optim.value_T34 <- opt$value
         T34$optim.counts_T34 <- opt$counts
         z$T34 <- T34

         #message('A=',A," K=",K," H=",H,"\n")
         if(! are.paraep4.valid(z)) {
            message <- "One or more parameters are not valid: One of the A methods."
            warning(message)
            z$ifailtext <- message
            z$ifail <- 3
         }
      #print(para)
      return(z)
}

# ifail3 is a parameter validity failure
# ifail2 is a general attempt to have a singular failure by sometype of eps outside of optim
# ifail1 is a failure by optim
