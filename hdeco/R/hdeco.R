"hdeco" <-
function (BE1=demoimage1, BE2=demoimage2, MICIKE=decomppath, MSK=FALSE, MASK=NULL, fnev="", AutoDecoPath=FALSE, JPG=FALSE, zsir=FALSE, MULTIX=FALSE, RECODEX=FALSE, NXRECODES=1, LUTX=NULL, RECODEZ=FALSE, NZRECODES=1, LUTZ=NULL, HISTOGRAM=FALSE, Z1DROP=FALSE, OMITX1=FALSE, PS=FALSE) {

  #############################################################
  # 
  # TITLE:  hdeco()
  # AUTHOR: TARMO REMMEL (BASED ON SANDOR KABOS)
  # DATE:   6 FEB, 2004
  # CALLS:  numerous...
  # NEEDS:  INPUT MATRIX IMAGE OBJECT(S), MASK IS OPTIONAL, RECODE LUT(S) ARE OPTIONAL, DECOMPOSITION PATH
  #
  # NOTES:      JPG = T|F      WHETHER IMAGES WILL BE SAVED AS JPEG
  #             zsir = T|F      WHETHER COLOUR DECOMPOSITION WILL BE CONDUCTED
  #             fnev = STRING   THE OUTPUT BASE FILENAME (WITHOUT EXTENSION)
  #             MICIKE = OBJ    MICIKE HAS PRIORITY OVER AutoDecoPath()
  #             MASK = OBJ      A MATRIX OF {0,1}, FOR NO MASK SET TO NULL
  #             BE1 = OBJ       INPUT IMAGE 1 (MUST EXIST)
  #             BE2 = OBJ       INPUT IMAGE 2 (OPTIONAL)
  #             MULTIX = T|F    BOOLEAN FLAG FOR INPUTTING OBJECT OF MULTIPLE IMAGE NAMES
  #                             FOR IMPORT - IT OVERRIDES BE1 AND BE2 WHICH SHOULD BE NULL
  #             XOBJ = OBJ      IF MULTIX = T, THEN XOBJ CONTAINS A VECTOR OF IMAGES TO
  #                             IMPORT FOR PROCESSING.  IT SHOULD LOOK LIKE THE FOLLOWING:
  #                             "1286806.e.60.20.txt" "1286776.e.60.20.txt" "1355636.e.60.20.txt"
  #             MSK = T|F       SHOULD THE AUTOMASK BE APPLIED?
  #				HISTOGRAM T|F	WHETHER HISTOGRAMS FOR ALL INPUT IMAGES WILL BE DRAWN ON A SECOND GRAPHICS DEVICE
  #				Z1DROP = T|F	IF RECODING Z VARIABLES... SHOULD THE ORIGINAL Z1 BE DROPPED? (SAVES MEMORY IF TRUE)
  #				OMITX1 = T|F	IF KNOWING THE UNIQUE IMAGE IDENTIFIERS IS NOT NECESSARY, THEY CAN BE OMITTED IN THE
  #								QKEP MATRIX CONSTRUCTION STAGE TO REDUCE MEMORY USAGE.  E.G., 324 NFI PLOTS ENTERING
  #                             INTO QKEP WOULD REQUIRE AN ADDITIONAL 324 DIMENSIONS TO THE QKEP ARRAY... HOWEVER, 
  #								WE'RE ALMOST ALWAYS MORE CONCERNED ABOUT THE X2, X3, ... VARIABLES, SO DROPPING X1
  #								MAKES LIFE SIMPLER.  IN REALITY, THE X1 VARIABLE STAYS IN THE MODEL BUT IS FILLED
  #								WITH 1s.  IT CAN BE ELIMINATED FROM ANALYSIS BY FILLING THE APPROPRIATE VFONAL COLUMN
  #								WITH -1s.
  #
  #             ENTROPY DECOMPOSITION WITH AUTOMATIC DECOMPOSITION PATH DEFINITION IF AutoDecoPath=TRUE, ELSE 
  #             MATRIX-LIKE DEFINITION OF DECOMPOSITION PATH IN MICIKE, OR, IF MICKIE=NULL - AN INTERACTIVE DEFINITION
  #
  #
  # OUTPUTS:    .QND            DIMENSION NAMES FOR .QKEP
  #             .QKEP           MULTIDIMENSIONAL ARRAY FOR HOLDING IMAGES
  #             .N              THE NUMBER OF PIXELS (ROWS IN KI MATRIX)
  #             .VFONAL         THE DECOMPOSITION PATH
  #				.NHIPO	NULL HYPOTHESES
  #				.AHIPO	ALTERNATE HYPOTHESES
  #             .MASKTITLE		TITLE OF THE GENERAL MASK USED
  #				.CIM	TITLES OF IMAGES PROCESSED
  #             
  #############################################################
  
  # DISPLAY STARTING MESSAGE AND PREPARE GRAPHICAL ENVIRONMENT
  cat("\nStarting Entropy Decomposition...")
  par(mfrow=c(3,3), pty="s")

  # STORE NUMBER OF IMAGES TO PROCESS
  NUMIMAGES <- 0
  
  # STORE GENERAL MASK TITLE IN .MASKTITLE IF GIVEN, OTHERWISE, SET TO NULL
  if(is.null(MASK)) {
    MASKTITLE <- NULL
    assign(".MASKTITLE", MASKTITLE, pos=1)
  }
  else {
    assign(".MASKTITLE", cim(MASK), pos=1)
  }
  
  # CREATE JPEG IF REQUIRED
  if(JPG){
    jpfnev <- ujfnev(fnev)
    # jpfnev <- paste(sep="",substring(jpfnev,1,nchar(jpfnev)-4),"a\%01d.jpeg")
    jpfnev <- paste(sep="",substring(jpfnev,1,nchar(jpfnev)-4),".jpeg")
    cat(jpfnev)
    jpeg(filename=jpfnev)
  }

  
  # ------------------------------------------------------
  # DECIDE WHICH INPUT AND .QKEP GENERATOR TO USE
  if(MULTIX) {
    cat("\nProcessing multiple images from vector object.")
    bendManyX4(INOBJ=XOBJ, GENERALMASK=MASK, AUTOMASK=MSK, XRECODE=RECODEX, NUMXRECODES=NXRECODES, XLUT=LUTX, ZRECODE=RECODEZ, NUMZRECODES=NZRECODES, ZLUT=LUTZ, DROP=Z1DROP, X1OMIT=OMITX1)
    NUMIMAGES <- length(XOBJ)
  }
  else {      # DO NOT READ MULTIPLE IMAGES FROM OBJECT
    if(is.null(BE2)) {
      cat("\nProcessing 1 image.")
      benc4(BE=BE1, GENERALMASK=MASK, AUTOMASK=MSK, XRECODE=RECODEX, NUMXRECODES=NXRECODES, XLUT=LUTX, ZRECODE=RECODEZ, NUMZRECODES=NZRECODES, ZLUT=LUTZ, DROP=Z1DROP, X1OMIT=OMITX1)
      NUMIMAGES <- 1
    }
    else {
      cat("\nProcessing 2 images.")
      bend4(BE1=BE1, BE2=BE2, GENERALMASK=MASK, AUTOMASK=MSK, XRECODE=RECODEX, NUMXRECODES=NXRECODES, XLUT=LUTX, ZRECODE=RECODEZ, NUMZRECODES=NZRECODES, ZLUT=LUTZ, DROP=Z1DROP, X1OMIT=OMITX1)
      NUMIMAGES <- 2
    }
  }
  # ------------------------------------------------------
  
  
  # DECIDE WHICH DECOMPOSITION PATH TO USE
  if(!is.null(MICIKE)) {
    cat("\nUsing custom decomposition path.")
    AutoDecoPath <- F
  }
  if(AutoDecoPath) {
    cat("\nUsing automatic decomposition path (Kabos).")
    MICIKE <- autodecopath()
  }
  vfonal(MICIKE)
  
     
  # PERFORM DECOMPOSITION
  cat("\nPerforming decomposition...\n\n")
  teszt(MICIKE,fnev,zsir)

  if(!MULTIX) {
    cat("\nGenerating graphical outputs...") 
    # GENERATE GRAPHICAL OUTPUTS
    if(is.null(BE2)){
      par(pty="s")
      if(!is.null(MASK)) {
        BE1 <- BE1 * MASK
      }
      imaks(BE1)
      title(sub="Masked Image")
      hprofil.rajz()
      par(pty="s")
      condhprofil.rajz()
      unchprofil.rajz()
    }
    else {
      if(!is.null(MASK)){ 
        BE1 <- BE1 * MASK
        BE2 <- BE2 * MASK
      }
      if(PS) {
        postscript(file=paste(substring(fnev,1,8),".eps", sep=""), onefile=T, paper="letter", horizontal=F)
      }
      par(pty="s", mfrow=c(2,2))
      fixedcolimage(BE1)
      title(sub="Image 1")
      fixedcolimage(BE2)
      title(sub="Image 2")
      condhprofil.rajz()
      unchprofil.rajz()
      if(PS) {
        dev.off()      
      }
    }
  }
  else {
    par(pty="s", mfrow=c(2,2))
    condhprofil.rajz()
    unchprofil.rajz()
  }
  
  
  # DRAW HISTOGRAMS FOR INPUT IMAGES IF DESIRED
  if(HISTOGRAM) {
  
    cat("\nGenerating histogram(s).")
  
    windows()
    par(mfrow=c(ceiling(NUMIMAGES / 2),2), pty="s")
        
    if(MULTIX) {
      for(rep in 1:NUMIMAGES) {
        hist(get(XOBJ[rep]), breaks=c(0:max(get(XOBJ[rep]))), xlab="Colour", right=TRUE, col=T, density=25, main=paste("Histogram for ", as.character(XOBJ[rep]), sep=""))
      }
    }
    else {
      if(!is.null(BE1)){
        hist(BE1, breaks=c(0:max(BE1)), xlab="Colour", right=TRUE, col=T, density=25)
      }
    
      if(!is.null(BE2)){
        hist(BE2, breaks=c(0:max(BE2)), xlab="Colour", right=TRUE, col=T, density=25)
      }
    }
    
  }
  
   
  # CLOSE GRAPHICS DEVICE IF JPEG OPTION IS TRUE
  if(JPG) {
    dev.off() 
  }

  return(cat("\nOK\n"))
}

