"condprofil.rajz" <-
function (BE=.HPROFIL) {

  #############################################################
  # 
  # TITLE:  unchprofil.rajz()
  # AUTHOR: TARMO REMMEL, (BASED ON SANDOR KABOS)
  # DATE:   16 JUNE, 2004
  # CALLS:  
  # CALLED BY:	t.hdeco() SERIES
  # NEEDS:  	.HPROFIL OBJECT
  # NOTES:  	GENERATES PLOT OF CONDITIONAL ENTROPY VERSUS DECOMPOSITION STEP
  #
  #############################################################

  # STORE THE NUMBER OF ROWS IN .HPROFIL
  NR <- dim(BE)[1]

  # GENERATE X DATA VECTOR OF INTEGERS FROM 1 TO NUMBER OF ROWS (NUMBER OF STEPS)
  X <- 1:NR

  # CREATE CONDITOINAL ENTROPY VECTOR
  # H(A|B) = H(A,B) - H(B)
  Y <- BE[,1] - BE[,3]

  # CREATE PLOT IF UNC COLUMN OF .HPROFIL DOES NOT CONTAIN ANY NaN ENTRIES
  if(TRUE %in% is.na(.HPROFIL[,5])) {
    cat("\nNaN in uncertainty results - Skipping uncertainty graph construction.")
  }
  else {
    sigplot(mat=.HPROFIL, column=4, sigcol=7, tit="H Decomposition", xtit="Step", ytit="Uncertainty Coefficient (%)", override=F, lowy=0, highy=1.0)
    #plot(x=X, y=Y, type="l", xaxt="n", xlab="Step", ylab="Entropy (NULL | ALT)", main="Conditional Entropy")
    #axis(1, at=X, labels=as.character(X))
  }
}

