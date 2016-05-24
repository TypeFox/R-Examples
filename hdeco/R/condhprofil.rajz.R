"condhprofil.rajz" <-
function (BE=.HPROFIL) {

  #############################################################
  # 
  # TITLE:  condhprofil.rajz()
  # AUTHOR: SANDOR KABOS
  # DATE:   20 AUG, 2003
  # CALLS:  
  # NEEDS:  
  # NOTES:  GENERATE PLOT OF MUTUAL INFORMATION VERSUS
  #			DECOMPOSITION STEP
  #
  #############################################################

  # THE FOLLOWING LINE IS A MANUAL FLAG TO TURN CONDITIONAL ENTROPY GRAPHING ON|OFF
  BOTH <- F
  
  
  # COMPUTE THE NUMBER OF ROWS IN .HPROFIL
  NR <- dim(BE)[1]
  
  # GENERATE THE X DATA VECTOR OF INTEGERS 1 TO NUMBER OF ROWS (STEPS)
  X <- 1:NR
  
  # EXTRACT THE Y DATA VECTOR OF MUTUAL INFORMATION VALUES
  YA <- BE[,4]
  
  # CREATE THE Y DATA VECTOR OF CONDITIONAL ENTROPIES
  # H(A|B) = H(A,B) - H(B)
  YB <- BE[,1] - BE[,3]
    
  # CREATE BLANK PLOT ENVIRONMENT  
  if(BOTH) {
    plot(x=X, y=YA, type="n", ylim=c(0, max(YA, YB)), xlab="Step", ylab="Entropy", main="H Decomposition")
  }
  else {
    par(pty="s")
    sigplot(mat=.HPROFIL, column=4, sigcol=7, tit="H Decomposition", xtit="Step", ytit="Mutual Information", override=F, lowy=0, highy=1.0)
    #plot(x=X, y=YA, type="n", xlab="Step", ylab="Mutual Information", main="H Decomposition")
  }
  
  # ADD MUTU AND CONDH PLOT LINES
  lines(x=X, YA, lty=1)
  if(BOTH) {
    lines(x=X, YB, lty=2)
    
    # GENERATE LEDGEND - x,y LOCATION IN GRAPH UNITS
    legend(NR - 3, max(YA, YB) * 0.75, legend = c("MUTU", "NULL | ALT "), lty=c(1,2), cex=0.9)
  }
    
}

