wi <-
function(BE=ClassPatternData$demoimage1, CONTROL=TRUE, PARAM1=ClassPatternData$WIBI, PARAM2=ClassPatternData$WITRUNC, solo=FALSE, SIZE=6) {

  #--------------------------------------------------------------
  # 
  # TITLE:     wi()
  # AUTHOR:    FERKO CSILLAG AND SANDOR KABOS, MODIFIED: TARMO REMMEL
  # DATE:      6 MAY 2013
  # CALLS:     NA
  # CALLED BY: wtest.loop()
  # NEEDS:     NA
  # NOTES:     VARIOUS HELPER CODE BITS THAT WORK WITH THE WHITTLE
  #            ESTIMATOR AND BIAS CORRECTOR
  #         
  #            TO RUN ON AN IMAGE WITHOUT THE LOOPING WHITTLE MATERIAL:
  #            rho <- wi(BE=demoimage2, solo=TRUE, SIZE=6)
  #
  #            OTHERWISE, THIS FUNCTION IS CALLED BY wtest.run(), 
  #            WHICH MAY BE CALLED BY build.lut()
  #
  #--------------------------------------------------------------

  wibi <- function(N=100) {
    KI <- rep(2, N)
    for(n in 2:N) {
      szor <- (2 * (2 * n - 1))/n
      KI[n] <- KI[n - 1] * szor
    }
    KI <- KI^2
    KI <- KI/c(1:N)
    return(list(KI, N))
  }

  wiga <- function(gamma=0.1) {
    g2 <- gamma^2
    KI <- array(dim = ClassPatternData$WITRUNC)
    KI[1] <- g2
    for(i in 2:ClassPatternData$WITRUNC) {
      KI[i] <- KI[i - 1] * g2
    }
    KI <- t(KI) %*% ClassPatternData$WIBI
    KI <- KI/2
    KI <- KI + log(CC$C0 - CC$C1 * gamma)
    return(KI)
  }

  wicc <- function(BE) {
    NR <- dim(BE)[1]
    NC <- dim(BE)[2]
    BE <- BE - mean(BE)
    BU <- c((BE[c(NR, (1:(NR - 1))),  ]))
    BF <- c((BE[, c((2:NC), 1)]))
    FUGG <- t(as.vector(BE)) %*% as.vector(BU)
    VIZSZ <- t(as.vector(BE)) %*% as.vector(BF)
    ITT <- t(as.vector(BE)) %*% as.vector(BE)
    return(list(C0 = ITT, C1 = 2 * (FUGG + VIZSZ)))
  }

  if(solo) {
    setuparray <- wibi(SIZE)
    assign("WIBI", setuparray[[1]], envir=ClassPatternData)
    assign("WITRUNC", setuparray[[2]], envir=ClassPatternData)
  }
  else {
    assign("WIBI", PARAM1, envir=ClassPatternData)
    assign("WITRUNC", PARAM2, envir=ClassPatternData)
  }

  CC <- wicc(BE)
  CO <- unlist(CC)[1]
  C1 <- unlist(CC)[2]
  KI <- nlminb(0.1, wiga, lower = 0, upper = 0.25)
  if(CONTROL) {
    write(c("rho , tau2 : ", round(KI$par, 4), round(CC$C0 - CC$C1 * KI$par, 4)))
  }
  GARB <- list(value = KI$par)

  return(GARB$value)
  
}
