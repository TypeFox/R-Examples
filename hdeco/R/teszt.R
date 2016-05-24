"teszt" <-
function (MICS=NULL,fnev="hh",zsir) {

  #############################################################
  # 
  # TITLE:  teszt()
  # AUTHOR: SANDOR KABOS, MODIFIED BY TARMO REMMEL
  # DATE:   27 NOV, 2003
  # CALLS:  ujfnev(), szabfok()
  # NEEDS:  
  # NOTES:  PERFORMS THE ENTROPY DECOMPOSITION AND WRITES LOG FILE
  #
  #############################################################

  # STORE THE DIMENSIONS OF .QKEP
  DIM <- dim(.QKEP)

  # STORE THE NUMBER OF ALTERNATE HYPOTHESES
  HHHOSSZ <- length(.AHIPO)

  # STORE JUST THE LETTER CODES FOR DIMENSIONS WITHOUT THE NUMERICAL SUBSCRIPTS
  DIMNAM1 <- substring(names(.QND), 1, 1)

  # SET THE NUMBER OF IMAGES
  EGYKE <- 1

  # IF THE FIRST DIMENSION OF .QKEP IS > 1, WE COMPARE MORE THAN ONE IMAGE (CURRENTLY MAX TWO)
  # CONSIDER SETTING EGYKE <- max(DIM[1])
  if(DIM[1] > 1) {
    EGYKE=2
  }

  # PREPARE MATRIX OF OUTPUT PHRASES
  NO <- matrix(c("(heterogeneous area)", "(homogeneous area)", "(different images)", "(similar images)"), nrow=2, byrow=TRUE)

  # GENERATE A UNIQUE OUTPUT FILENAME
  fnev <- ujfnev(fnev)
  if(fnev == "") {
    tempfnev <- "- NA -"
  }
  else {
    tempfnev <- substring(fnev, 1, nchar(fnev) - 4)
  }
  
    
  # WRITE LOG HEADER INFORMATION ABOUT DATE, FILENAME, AND MASK USED
  #write(c("Entropy Decomposition,", date(), ", file name:", substring(fnev,1,nchar(fnev)-4)), ncol=55, file=fnev)
  write(c("Entropy Decomposition,",date(),", File name:",tempfnev), ncol=55, file=fnev)
  if(is.null(.MASKTITLE)) {
    .MASKTITLE <- "No general mask applied."
  }
  write(.MASKTITLE, ncol=55, file=fnev, append=TRUE)

  # WRITE LOG HEADER INFORMATION ABOUT NUMBER OF IMAGES
  if(EGYKE==1) {
    write("One image", ncol=55, file=fnev, append=TRUE)
  }
  else {
    XEK <- "X" == DIMNAM1
    #NIX<-prod(DIM[XEK])  # KABOS VERSION (PRODUCT)
    NIX <- max(DIM[XEK])
    write(c(NIX, "images, coded as X variable(s)"), ncol=55, file=fnev, append=TRUE)
  }

  # COMPUTE AND STORE BINARY VECTOR INDICATING TRUE FOR Y VARIABLE DIMENSIONS
  YEK <- "Y"==DIMNAM1
  # STORE PRODUCT OF Y DIMENSIONS
  NIY <- prod(DIM[YEK])
  # STORE THE NUMBER OF Y VARIABLES
  NLY <- sum(YEK)
  # STORE OUTPUT PHRASES DEPENDING ON NUMBER OF IMAGES 
  NEVELO <- c("The","Each")
  # WRITE TO LOG ABOUT THE NUMBER OF IMAGE PIXELS AND Y vARIABLES  
  write(c(NEVELO[EGYKE],"image consists of",NIY,"pixels in a",NLY,"-level pyramid, coded as Y variable(s)"), ncol=55, file=fnev, append=TRUE)

  # COMPUTE AND STORE BINARY VECTOR INDICATING TRUE FOR Z VARIABLE DIMENSIONS
  ZEK <- "Z" == DIMNAM1
  # STORE PRODUCT OF Z DIMENSIONS  
  NIZ <- prod(DIM[ZEK])
  
  # WRITE TO LOG ABOUT COLOUR CODES
  # ORIGINAL LINE IS COMMENTED OUT - FOLLOWED BY TARMO'S MODIFICATION FOR NUMBER OF COLOURS
  #write(c(NIZ,"colors, coded as Z variable(s)\n"), ncol=55, file=fnev, append=TRUE)
  write(c(.COLOURS,"colors, coded as Z variable(s)"), ncol=55, file=fnev, append=TRUE)
  
  
  # WRITE TO LOG ABOUT THE TITLE OF IMAGE(S) AND DECOMPOSITION PATH
  if(is.null(.CIM)) {
    .CIM <- "No image title(s) given."
  }
  write(.CIM, ncol=55, file=fnev, append=TRUE)
  write(cim(.VFONAL), ncol=55, file=fnev, append=TRUE)

  # CONSTRUCT THE ENTROPY PROFILE SUMMARY MATRIX
  #HPROFIL <- matrix(nrow=HHHOSSZ,ncol=5,dimnames=list(NULL,c("HALAPF","HNULL","HALT","MUTU","UNC")))
  # NEXT LINE IS TARMO'S UPDATE TO ADD G^2 AND SIGNIFICANCE VALUES TO .HPROFIL
  HPROFIL <- matrix(nrow=HHHOSSZ, ncol=9, dimnames=list(NULL,c("HALAPF","HNULL","HALT","MUTU","UNC","G2","SIG","SING-MULT","DESC")))

  # FOR EACH RECORD IN THE ALTERNATIVE HYPOTHESIS MATRIX
  for(ik in 1:HHHOSSZ) { 
    KELLFELT <- !is.null(.KIVALO[[ik]])
    KI <- tesz(ik)
    MUTU <- KI$HHIPO + KI$HAHIPO - KI$HALAPF
    UNC <- 100* MUTU / KI$HAHIPO
    TESZT <- MUTU * 2 * .N
    if(abs(TESZT) < 1e-7) {     # IF tesz() INDICATES ERRORS - THEY CAUSE CRASHES HERE
      TESZT <- 0
    }
    DF <- szabfok(KI$ALAPF) - szabfok(.NHIPO[[ik]]) - szabfok(.AHIPO[[ik]]) 

    # COMPUTE SIGNIFICANCE DEPENDING ON DEGREES OF FREEDOM    
    if(DF <= 0) {
      SIG <- 1
    }
    else {
      SIG <- pchisq(TESZT, DF, , FALSE)
    }

    # SET NINI TO GOVERN PHRASE OUTPUT TO LOG LATER ON DEPENDING ON SIGNIFICANCE VALUE    
    if(SIG < 0.05) {
      NINI <- 1
    }
    else {
      NINI <- 2
    }
    
    # UPDATE THE ENTROPY PROFILE MATRIX FOR THE GIVEN STEP
    #HPROFIL[ik,] <- c(KI$HALAPF,KI$HHIPO,KI$HAHIPO,MUTU,UNC)
    # NEXT LINE IS TARMO'S UPDATE
    HPROFIL[ik,] <- c(KI$HALAPF, KI$HHIPO, KI$HAHIPO, MUTU, UNC, TESZT, SIG, EGYKE, NINI)
    
    # SET STEP NUMBER FOR LOG OUTPUT
    if(EGYKE == 1) {
      ikstep <- ik
    }
    else {
      ikstep <- ik - 1
    }

    # WRITE DECOMPOSITION RESULTS TO LOG
    write(c("\nStep", ikstep, "\nBaseH:", names(KI$ALAPF), ",  H=", signif(KI$HALAPF,5), ",       DF=", szabfok(KI$ALAPF)), ncol=55, file=fnev, append=TRUE)
    write(c("NullH:", names(.NHIPO[[ik]]), ",  H=", signif(KI$HHIPO,5), ",       DF=", szabfok(.NHIPO[[ik]])), ncol=55, file=fnev, append=TRUE)
    write(c("Alt H:", names(.AHIPO[[ik]]),",  H=",signif(KI$HAHIPO,5),",       DF=",szabfok(.AHIPO[[ik]])), ncol=55, file=fnev, append=TRUE)
    write(c("Mutual info (", names(.AHIPO[[ik]]), "between", names(.NHIPO[[ik]]), ") =", signif(MUTU,5), ",       DF=",    DF), ncol=55, file=fnev, append=TRUE)
    write(c("Uncertainty coeff (", names(.AHIPO[[ik]]), "given", names(.NHIPO[[ik]]), ") =", paste(sep="",signif(UNC,4),"%")), ncol=55, file=fnev, append=TRUE)
    write(c("Test statistics (", names(.AHIPO[[ik]]), "joint) =", signif(TESZT,5), ",   DF=", DF, ",  signif =", signif(SIG,5), NO[EGYKE,NINI]), ncol=55, file=fnev, append=TRUE)
    
    if(KELLFELT) {
      COMUTU <- KI$HHIPO + KI$HBEVALO.AHIPO - KI$HBEVALO.ALAPF
      COMUTU <- MUTU - COMUTU
      COTESZT <- COMUTU*2*.N
      CODF <- szabfok(KI$BEVALO.ALAPF) - szabfok(.NHIPO[[ik]]) - szabfok(KI$BEVALO.AHIPO) 
      CODF <- DF - CODF
      COSIG <- pchisq(COTESZT,CODF,,F)
      if(COSIG < 0.05) {
        CONINI <- 1
      }
      else {
        CONINI <- 2
      }
      # MAKE LOG ENTRY
      write(c("Conditional H:  (",names(.KIVALO[[ik]]),"conditioned on",names(KI$BEVALO.AHIPO),") =",signif(COMUTU,5)), ncol=55, file=fnev, append=TRUE)
    }
  }
  
  # WRITE THE .HPROFIL TABLE TO THE EXPORT FILE
  write("\n\n.HPROFIL Table with headers", file=fnev, append=TRUE)
  write.table(HPROFIL, file=fnev, sep=",", append=TRUE)
  
  assign(".HPROFIL", HPROFIL, pos=1)
  
  # PERFORM COLOUR DECOMPOSITION IF zsir IS TRUE
  if(zsir) {
    cat("\nPerforming colour decomposition (zsir)...")
    zsir(, fnev, EGYKE)
  }
  
  return(cat("\nDecomposition complete."))
}

