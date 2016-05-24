"pwm2lmom" <-
function(pwm) {

    if(is.list(pwm)) {
      # are the betas separate and we then expect all five?
      # or are the betas available as a vector embedded in the list?
      if(! is.null(pwm$BETA0)) {
        z <- list(L1   = NULL,
                  L2   = NULL,
                  TAU3 = NULL,
                  TAU4 = NULL,
                  TAU5 = NULL,
                  LCV  = NULL,
                  L3   = NULL,
                  L4   = NULL,
                  L5   = NULL
                 )
        z$L1 <- pwm$BETA0
        z$L2 <- 2*pwm$BETA1 - pwm$BETA0
        z$L3 <- 6*pwm$BETA2 - 6*pwm$BETA1 + pwm$BETA0
        z$L4 <- 20*pwm$BETA3 - 30*pwm$BETA2 + 12*pwm$BETA1 - pwm$BETA0
        z$L5 <- 70*pwm$BETA4 - 140*pwm$BETA3 + 90*pwm$BETA2 - 20*pwm$BETA1 + pwm$BETA0
        z$LCV <- z$L2/z$L1
        z$TAU3 <- z$L3/z$L2
        z$TAU4 <- z$L4/z$L2
        z$TAU5 <- z$L5/z$L2
        return(z)
      }
      # betas are in a list element, reset them to the pwm for the
      # remainder of the function to operate
      if(! is.null(pwm$betas)) {
        pwm <- pwm$betas
      }
      else {
        warning("ambiguous call, do not find Betas for processing")
        return(NULL)
      }
    }

    nmom <- length(pwm)

    L <- vector(mode = "numeric",length=nmom)
    R <- vector(mode = "numeric",length=nmom)
      
    for(i in seq(1,nmom)) {
      r <- i - 1
      sum <- 0
      for(k in seq(0,r)) {
        weight <- (-1)^(r-k)*choose(r,k)*choose(r+k,k)
        sum <- sum + weight*pwm[k+1]
      }
      L[i] <- sum
    }
    if(nmom >= 2) {
      R[2] <- L[2]/L[1]
    }
    if(nmom >= 3) {
      for(r in seq(3,nmom)) {
        R[r] <- L[r]/L[2]
      }
    }
    R[1] <- NA

    z <- list(lambdas=L,
              ratios=R,
              source="pwm2lmom")
    return(z)
}

