wtest.run <-
function(LEVEL=6, REPSIM=20, RHO=0.2499999, CPROP=0.5, RAJZ=T, CIM="CIM") {

  #--------------------------------------------------------------
  # 
  # TITLE:     wtest.run()
  # AUTHOR:    FERKO CSILLAG, MODIFIED BY TARMO REMMEL
  # DATE:      15 AUGUST 2003, 16 JULY 2013
  # CALLS:     wibi(), CARsimu(), wi() --> wicc()
  # CALLED BY: NA
  # NEEDS:     NA
  # NOTES:     THIS PROGRAM USES THE FFT-METHOD TO SIMULATE IMAGES
  #            WITH A GIVEN LEVEL OF RHO.  THEN THE (ISOTROPIC,
  #            FIRST-NEIGHBOUR AUTOCORRELATION PARAMETER IS
  #            ESTIMATED TWICE: (1) ON THE ACTUAL CONTINUOUS
  #            SIMULATED IMAGE, AND (2) ON A "CUT" BINARY IMAGE,
  #            WHERE THE PROPORTION OF WHITE (0) IS DEFINED BY
  #            CPROP USING P. WHITTLE'S (1954!) METHOD (FOR
  #            REGULAR GRIDS) [ON STATIONARY PROCESSES IN THE PLANE.
  #            BIOMETRIKA 41:434-449.] AND THE RESULTS ARE SUMMARIZED.
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

  RESULT <- rep(0, 2 * REPSIM)
  dim(RESULT) <- c(2, REPSIM)

  # SETUP ARRAYS FOR THE WHITTLE ESTIMATION
  setuparray <- wibi(LEVEL)
  assign("WIBI", setuparray[[1]], envir=ClassPatternData)
  assign("WITRUNC", setuparray[[2]], envir=ClassPatternData)

  for(lup in 1:REPSIM) {
    W <- CARsimu(rho = RHO, rajz = F)
    RESULT[1, lup] <- wi(BE = W, CONTROL = RAJZ, PARAM1 = ClassPatternData$WIBI, PARAM2 = ClassPatternData$WITRUNC, solo=FALSE, SIZE=LEVEL)
    TEMP <- quantile(W, CPROP)
    GARB <- W > TEMP[1]
    GARB <- factor(GARB)
    GARB <- as.numeric(GARB)
    W.0 <- GARB
    dim(W.0) <- c(2^LEVEL, 2^LEVEL)
    if(RAJZ) {
      par(mfrow = c(1, 2))
      par(pty = "s")
      image(W)
      image(W.0)
      cat("\r... PRESS ENTER TO CONTINUE (1 TO LET IT RUN)...")
      ans <- readline()
      if(ans == 1)
      RAJZ <- F
    }
    RESULT[2, lup] <- wi(BE = W.0, CONTROL = RAJZ, PARAM1 = ClassPatternData$WIBI, PARAM2 = ClassPatternData$WITRUNC, solo=FALSE, SIZE=LEVEL)
    cat("\r                          ", lup, "ITERATION OUT OF:", REPSIM)
  }
  RESULT <- unlist(RESULT)
  dim(RESULT) <- c(2, REPSIM)

  # PLOTTING
  par(mfrow = c(1, 1))
  par(pty = "s")
  if(CIM != "") {
    boxplot(RESULT[1,  ], RESULT[2,  ], names = c("ORIGINAL-CONTINUOUS", "BINARY"), ylim = c(0, 0.25))
    mastertitle <- paste("RHO | PROP : ", as.character(RHO), " | ", as.character(CPROP), sep = " ")
    title(mastertitle)
  }
  
  return(RESULT)
  	
}
