"unchprofil.rajz" <-
function (BE=.HPROFIL) {

  #############################################################
  # 
  # TITLE:  unchprofil.rajz()
  # AUTHOR: SANDOR KABOS, MODIFIED BY TARMO REMMEL
  # DATE:   27 NOV, 2003
  # CALLS:  
  # CALLED BY:	t.hdeco() SERIES
  # NEEDS:  	.HPROFIL OBJECT
  # NOTES:  	GENERATES PLOT OF UNCERTAINTY VERSUS DECOMPOSITION STEP
  #
  #############################################################

  # STORE THE NUMBER OF ROWS IN .HPROFIL
  NR <- dim(BE)[1]

  # GENERATE X DATA VECTOR OF INTEGERS FROM 1 TO NUMBER OF ROWS (NUMBER OF STEPS)
  X <- 1:NR

  # EXTRACT Y DATA VECTOR OF UNCERTAINTY (%)
  Y <- BE[,5]

  # CREATE PLOT IF UNC COLUMN OF .HPROFIL DOES NOT CONTAIN ANY NaN ENTRIES
  if(TRUE %in% is.na(.HPROFIL[,5])) {
    cat("\nNaN in uncertainty results - Skipping uncertainty graph construction.")
  }
  else {
    par(pty="s")
    sigplot(mat=.HPROFIL, column=5, sigcol=7, tit="Uncertainty", xtit="Step", ytit="Uncertainty Coefficient (%)", override=F, lowy=0, highy=1.0)
    #plot(x=X, y=Y, type="l", xaxt="n", xlab="Step", ylab="Uncertainty Coefficient (%)", main="Uncertainty")
    #axis(1, at=X, labels=as.character(X))
  }
}

