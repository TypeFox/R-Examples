findrow <-
function(autocorr = 0.2, VERBOSE=FALSE) {

  #--------------------------------------------------------------
  # 
  # TITLE:     findrow()
  # AUTHOR:    TARMO REMMEL
  # DATE:      16 JULY 2013
  # CALLS:     N/A
  # CALLED BY: singlemap()
  # NEEDS:     LOOKUP TALBE DIFF50 FOR BIAS CORRECTION AND
  #            A RHO VALUE TO LOOK UP
  # NOTES:     USED TO FIND THE CLOSEST VALUE IN THE DIFF50
  #            LOOKUP TABLE MATRIX FOR RHO BIAS CORRECTION
  #            USED IN CONJUNCTION WITH findcol()
  #--------------------------------------------------------------

  # RUN ONLY IF THE AUTOCORRELATION IS > 0 and < 0.2499999
  if(autocorr > 0 & autocorr < 0.2499999) {

    minrow <- max(which(rownames(ClassPatternData$DIFF50) < autocorr))
    maxrow <- min(which(rownames(ClassPatternData$DIFF50) > autocorr))

    val11 <- abs(as.numeric(rownames(ClassPatternData$DIFF50)[minrow]) - autocorr)
    val22 <- abs(as.numeric(rownames(ClassPatternData$DIFF50)[maxrow]) - autocorr)

    if(VERBOSE) {
      cat("minrow: ", minrow, "\n", sep="")    
      cat("maxrow: ", maxrow, "\n", sep="")
      cat("val12: ", val11, "\n", sep="")    
      cat("val22: ", val22, "\n", sep="")
    }

    if(val11 < val22) row <- as.numeric(minrow) else row <- as.numeric(maxrow)
  }
  else {
    # A SPECIAL CASE, INDICATING THAT THE CORRECTION IS NEGLIGIBLE AT THIS LEVEL
    row <- 99
  }
  
  return(row)

}
