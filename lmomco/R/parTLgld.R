"parTLgld" <-
function(lmom, verbose=FALSE, initkh=NULL, eps=1e-3,
         aux=c("tau5", "tau6"), checklmom=TRUE) {
    aux <- match.arg(aux)

    if(is.null(lmom$trim)) {
      warning("An L-moment object and not a TL-moment object appears to have been passed")
      return()
    }
    if(length(lmom$trim) == 1 && lmom$trim != 1) {
      warning("Attribute of TL-moments is not trim=1--can not complete parameter estimation")
      return()
    }

    LM1 <- lmom$lambdas[1]
    LM2 <- lmom$lambdas[2]

    T3 <- lmom$ratios[3]
    T4 <- lmom$ratios[4]
    T5 <- lmom$ratios[5]
    T6 <- lmom$ratios[6]

    # Four parameter distributions do not normally need the
    # fifth L-moment, but for the GLD as implemented here--we do
    # This error message is to help with this fact.
    if(aux == "tau5" && is.na(T5)) {
      warning("The fifth L-moment ratio TAU5 is undefined")
      return()
    }
    if(aux == "tau6" && is.na(T6)) {
      warning("The sixth L-moment ratio TAU6 is undefined")
      return()
    }

    # THE TRIMMED L-MOMENTS OF GLD
    estTLla1 <- function(La2,La3,La4) {
      La1 <- LM1 - 6*La2*(1/((La3+3)*(La3+2)) - 1/((La4+3)*(La4+2)))
    }

    estTLla2 <- function(LM2,La3,La4) {
      D1 <- (La3+2)*(La3+3)*(La3+4)
      D2 <- (La4+2)*(La4+3)*(La4+4)
      return((LM2/6)/( La3/D1 + La4/D2 ) )
    }

    estTLtau3 <- function(K,H) {
      D1 <- (K+5)*(K+4)*(K+3)*(K+2)
      D2 <- (H+5)*(H+4)*(H+3)*(H+2)
      G  <- K*(H+4)*(H+3)*(H+2) + H*(K+4)*(K+3)*(K+2)
      t3 <- (10/9)*(K*(K-1)*D2 - H*(H-1)*D1)/((K+5)*(H+5)*G)
      return(t3)
    }
    estTLtau4 <- function(K,H) {
      D1 <- (K+6)*(K+5)*(K+4)*(K+3)*(K+2)
      D2 <- (H+6)*(H+5)*(H+4)*(H+3)*(H+2)
      G  <- K*(H+4)*(H+3)*(H+2)+H*(K+4)*(K+3)*(K+2)
      t4 <- (5/4)*(K*(K-1)*(K-2)*D2 + H*(H-1)*(H-2)*D1) /
            ((K+6)*(H+6)*(K+5)*(H+5)*G)
      return(t4)
    }
    estTLtau5 <- function(K,H) {
      D1 <- (K+7)*(K+6)*(K+5)*(K+4)*(K+3)*(K+2)
      D2 <- (H+7)*(H+6)*(H+5)*(H+4)*(H+3)*(H+2)
      G  <- K*(H+4)*(H+3)*(H+2)+H*(K+4)*(K+3)*(K+2)
      t5 <- (7/5)*(K*(K-1)*(K-2)*(K-3)*D2 - H*(H-1)*(H-2)*(H-3)*D1)/
                ((K+7)*(H+7)*(K+6)*(H+6)*(K+5)*(H+5)*G)
      return(t5)
    }
    estTLtau6 <- function(K,H) {
      D1 <- (K+8)*(K+7)*(K+6)*(K+5)*(K+4)*(K+3)*(K+2)
      D2 <- (H+8)*(H+7)*(H+6)*(H+5)*(H+4)*(H+3)*(H+2)
      G  <- K*(H+4)*(H+3)*(H+2)+H*(K+4)*(K+3)*(K+2)
      t6 <- (14/9)*(K*(K-1)*(K-2)*(K-3)*(K-4)*D2 + H*(H-1)*(H-2)*(H-3)*(H-4)*D1) /
                   ((K+8)*(H+8)*(K+7)*(H+7)*(K+6)*(H+6)*(K+5)*(H+5)*G)
      return(t6)
    }
    # END OF TRIMMED L-MOMENTS OF GLD

    # ORDINARY TAU3, TAU4, AND TAU5 OF GLD GIVEN PARAMETERS
    esttau3 <- function(La3,La4) {
      N1 <- La3*(La3-1)*(La4+3)*(La4+2)*(La4+1)
      N2 <- La4*(La4-1)*(La3+3)*(La3+2)*(La3+1)
      D1 <- (La3+3)*(La4+3)
      D2 <- La3*(La4+1)*(La4+2) + La4*(La3+1)*(La3+2)
      t3 <- (N1 - N2)/(D1*D2)
      return(t3)
    }
    esttau4 <- function(La3,La4) {
      N1 <- La3*(La3-2)*(La3-1)*(La4+4)*(La4+3)*(La4+2)*(La4+1)
      N2 <- La4*(La4-2)*(La4-1)*(La3+4)*(La3+3)*(La3+2)*(La3+1)
      D1 <- (La3+4)*(La4+4)*(La3+3)*(La4+3)
      D2 <- La3*(La4+1)*(La4+2) + La4*(La3+1)*(La3+2)
      t4 <- (N1 + N2)/(D1*D2)
      return(t4)
    }
    esttau5 <- function(La3,La4) {
      N1 <- La3*(La3-3)*(La3-2)*(La3-1)*(La4+5)*(La4+4)*(La4+3)*(La4+2)*(La4+1)
      N2 <- La4*(La4-3)*(La4-2)*(La4-1)*(La3+5)*(La3+4)*(La3+3)*(La3+2)*(La3+1)
      D1 <- (La3+5)*(La4+5)*(La3+4)*(La4+4)*(La3+3)*(La4+3)
      D2 <- La3*(La4+1)*(La4+2) + La4*(La3+1)*(La3+2)
      t5 <- (N1 - N2)/(D1*D2)
      return(t5)
    }
    esttau6 <- function(La3,La4) {
      N1 <- La3*(La3-4)*(La3-3)*(La3-2)*(La3-1)*(La4+6)*(La4+5)*(La4+4)*(La4+3)*(La4+2)*(La4+1)
      N2 <- La4*(La4-4)*(La4-3)*(La4-2)*(La4-1)*(La3+6)*(La3+5)*(La3+4)*(La3+3)*(La3+2)*(La3+1)
      D1 <- (La3+6)*(La4+6)*(La3+5)*(La4+5)*(La3+4)*(La4+4)*(La3+3)*(La4+3)
      D2 <- La3*(La4+1)*(La4+2) + La4*(La3+1)*(La3+2)
      t6 <- (N1 + N2)/(D1*D2)
      return(t6)
    }
    # END ORDINARY L-MOMENTS

    # Define the objective function
    fn <- function(x) {
      La3 <- x[1]
      La4 <- x[2]
      t3  <- estTLtau3(La3,La4) # TRIMMED LMOMENTS USED FOR SOLUTION
      # BUT CONVERSION LATER ON TO ORDINARY USED TO TEST LMOMENT VALIDITY
      t4  <- estTLtau4(La3,La4)
      ss  <- ((T3-t3)^2 + (T4-t4)^2)
      return(ss)
    }

    # This function is a scaled down version of are.pargld.valid
    # to avoid the unnecessary overhead of computing the first
    # L-moment and building the standard parameter object.
    validgld <- function(La2,La3,La4) {
      if(is.na(La2)) return(FALSE)
      if(La2 == -Inf || La2 == Inf) return(FALSE)
      #if(verbose == TRUE) cat(c("validgld-Checking Theorem 1.3.3: ",La2, La3, La4,"\n"))
      # Test that second L-moment is suitable
      for(F in seq(0,1,by=0.00001)) {
        tmp <- La2*(La3*F^(La3-1) + La4*(1-F)^(La4-1))
        #cat(c(La2,La3,La4,tmp))
        if(tmp < 0) return(FALSE)
      }

      #if(verbose == TRUE) cat(c("validgld-Checking by region: ",La2, La3, La4,"\n"))
      # ratios define the curved lines in figure1.3-1 of K&D
      ratio6 <- -1/La3
      ratio5 <-  1/La4
      # See Theorem 1.3.33 of Karian and Dudewicz and figure1.3-1
      if(La3 <= -1 && La4 >=  1) {         # REGION 1
        return(FALSE) # ordinary L-moments do not exist in REGION 1
        return(TRUE)
      }
      else if(La3 >=  1 && La4 <= -1) {    # REGION 2
        return(FALSE) # ordinary L-moments do not exist in REGION 2
        return(TRUE)
      }
      else if(La3 >= 0 && La4 >= 0) {      # REGION 3
        return(TRUE)
      }
      else if(La3 <= 0 && La4 <= 0) {      # REGION 4
        return(TRUE)
      }
      else if((La4 < ratio6 && La4 >= -1) && La3 >= 1) {  # REGION 6
        return(TRUE)
      }
      else if(La4 > ratio5 && (La3 >= -1 && La3 <= 0)) { # REGION 5
        return(TRUE)
      }
      return(FALSE)
    }

    # Are the L-moments valid?  Not fully certain that this check
    # is needed if the GLD proves to be valid through the validgld
    # function above. Early testing of pargld() before the suitabiliy
    # of the second L-moment for GLD showed that the validlmom() test
    # was needed.  Consider this for safety at any rate.
    validlmom <- function(sol,attempt) {
      LM3 <- sol$par[1]
      LM4 <- sol$par[2]
      # TL-moment (t=1) based GLD only valid for > -2 parameters--Hosking email
      if(LM3 <= -2 || LM4 <= -2) return(FALSE)
      T3 <- esttau3(LM3,LM4)
      T4 <- esttau4(LM3,LM4)
      T5 <- esttau5(LM3,LM4)
      #if(verbose == TRUE) cat(c("validlmom(region,Tau3,Tau4,Tau5)?:\n",attempt,T3,T4,T5))
      if(abs(T3) > 1)  return(FALSE)
      if(T4 < (0.25*(5*T3^2 - 1)) || T4 > 1) return(FALSE)
      if(abs(T5) > 1)  return(FALSE)
      return(TRUE)
    }

   if(is.null(initkh)) {
     g <- 13
     M <- matrix(nrow = g, ncol = 2)
     M[1,]  <- c(-2.5,2.5)   # REGION 1
     M[2,]  <- c(2.5,-2.5)   # REGION 2
     M[3,]  <- c(2.5,2.5)    # REGION 3
     M[4,]  <- c(-.5,-.5)    # REGION 4
     M[5,]  <- c(-1.5,-1.5)  # REGION 4
     M[6,]  <- c(-2.5,-2.5)  # REGION 4
     M[7,]  <- c(-3.5,-3.5)  # REGION 4
     M[8,]  <- c(-4.5,-4.5)  # REGION 4
     M[9,]  <- c(-5.5,-5.5)  # REGION 4
     M[10,] <- c(-6.5,-6.5)  # REGION 4
     M[11,] <- c(-7.5,-7.5)  # REGION 4
     M[12,] <- c(-.5,2.5)    # REGION 5
     M[13,] <- c(2.5,-.5)    # REGION 6
   }
   else {
     M <- matrix(nrow = 1, ncol = 2)
     M[1,] <- initkh
   }

   num_eachs <- length(M[,1])

   each_attempt  <- vector(mode = 'numeric', length = num_eachs)
   each_initialK <- each_initialH <- each_error <- each_attempt
   each_interpretederror <- each_xi <- each_alpha <- each_attempt
   each_kappa <- each_h <- each_tdiff <- each_attempt
   each_valid_gld <- each_valid_lmom <- each_fine <-
       vector(mode = 'logical', length = num_eachs)

   for(i in seq(1,nrow(M))) {
     if(verbose) cat("  ",i,"(",num_eachs,")...", sep="")
     # Test for NaN from the KEK guess
     if(is.nan(M[i,1]) || is.nan(M[i,2])) next

     r       <- optim(M[i,],fn) # THE MAGIC IS HERE!!!!!!!!!!
     e       <- r$value # extract error
     K       <- r$par[1] # extract the two solutions
     H       <- r$par[2]
     TLdiff  <- ifelse(aux == "tau6",  estTLtau6(K,H) - T6,
                                       estTLtau5(K,H) - T5)
     # compute difference
     A       <- estTLla2(LM2,K,H)
     X       <- estTLla1(A,K,H)

     valid.lmom  <- ifelse(validlmom(r,attempt=i), TRUE, FALSE)
     valid.gld   <- ifelse(validgld(A,K,H),        TRUE, FALSE)
     good.enough <- ifelse(e < eps,                TRUE, FALSE)

     each_attempt[i]    <- i
     each_initialK[i]   <- M[i,1]
     each_initialH[i]   <- M[i,2]
     each_xi[i]         <- X
     each_alpha[i]      <- A
     each_kappa[i]      <- K
     each_h[i]          <- H
     each_tdiff[i]      <- TLdiff
     each_error[i]      <- e
     each_valid_gld[i]  <- valid.gld
     each_valid_lmom[i] <- valid.lmom
     each_fine[i]       <- good.enough
   }

   EACH <- data.frame(xi         = each_xi,
                      alpha      = each_alpha,
                      kappa      = each_kappa,
                      h          = each_h,
                      attempt    = each_attempt,
                      delTauHI   = each_tdiff,
                      error      = each_error,
                      initial_k  = each_initialK,
                      initial_h  = each_initialH,
                      valid.gld  = each_valid_gld,
                      valid.lmr  = each_valid_lmom,
                      lowerror   = each_fine)
   if(verbose) print(EACH)

   EACH <- EACH[EACH$valid.gld == TRUE,]
   EACH <- EACH[EACH$valid.lmr == TRUE,]
   EACH <- EACH[order(EACH$error),]

   BEST <- EACH[EACH$lowerror == TRUE,]
   BEST <- BEST[order(abs(BEST$delTauHI), BEST$error),]

   REST <- BEST[-c(1),]
   rownames(REST) <- NULL
   # Preparing final best guess . . .
   para <- vector(mode="numeric", length=4)
   names(para) <- c("xi","alpha","kappa","h")
   para[1]  <- BEST$xi[1]
   para[2]  <- BEST$alpha[1]
   para[3]  <- BEST$kappa[1]
   para[4]  <- BEST$h[1]
   taudiff  <- BEST$delTauHI[1]
   error    <- BEST$error[1]
   n <- length(REST[,1])
   if(n == 0) REST <- NULL
   return(list(type     = 'gld',
               para     = para,
               delTauHI = taudiff,
               error    = error,
               source   = "parTLgld",
               rest     = REST))
}
