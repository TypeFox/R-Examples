build.lut <-
function(LEVEL=6, REPSIM=5, RAJZ=FALSE, CIM="") {

  #--------------------------------------------------------------
  # 
  # TITLE:     build.lut()
  # AUTHOR:    FERKO CSILLAG, MODIFIED BY TARMO REMMEL
  # DATE:      15 AUGUST 2003, 11 MAY 2011
  # CALLS:     wtest.run()
  # CALLED BY: NA
  # NEEDS:     NA
  # NOTES:     USES wtest.run() TO RUN A LOOP ALONG AUTOCORRELATION
  #            (rho) AND PROPORTION AND THE RESULTS ARE SUMMARIZED
  #            AS A LOOKUP TABLE FOR BIAS ADJUSTMENT ON RHO
  #            ESTIMATION ON A BINARY MAP.
  #--------------------------------------------------------------

  DIFF <- rep(0, 110)
  dim(DIFF) <- c(10, 11)
  IX <- 0
  IY <- 0
  for(luprho in seq(0, 0.2499999, 0.0277777)) {
    RHO <- luprho
    IX <- IX + 1
    IY <- 0
    for(lupcprop in seq(0.1, 0.9, 0.1)) {
      CPROP <- lupcprop
      IY <- IY + 1
      for(lup in 1:REPSIM) {
        RESULTT <- wtest.run(REPSIM = REPSIM, LEVEL = LEVEL, RHO= RHO, CPROP = CPROP, RAJZ = RAJZ, CIM = CIM)
      }
      DIFF[IX, IY] <- median(RESULTT[1,  ]) - median(RESULTT[2,])
      }
  }

  # RETURN THE BIAS CORRECTION MATRIX TO THE USER
  return(DIFF)

}
