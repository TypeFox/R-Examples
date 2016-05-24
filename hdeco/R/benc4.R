"benc4" <-
function (BE=KEPP, GENERALMASK=NULL, AUTOMASK=TRUE, XRECODE=FALSE, NUMXRECODES=1, XLUT=NULL, ZRECODE=FALSE, NUMZRECODES=1, ZLUT=NULL, DROP=TRUE, X1OMIT=TRUE) {
  
  #############################################################
  # 
  # TITLE:  t.benc4()
  # AUTHOR: TARMO REMMEL
  # DATE:   18 JAN, 2004
  # CALLS:  benya(), kateg(), benci()
  # NEEDS:  BEMT, INPUT IMAGE MATRIX OBJECT, MASK IS OPTIONAL
  # NOTES:  READS A SINGLE IMAGE INTO A MULTIDIMENSIONAL ARRAY
  #         COLOUR CODES MUST BE 1,2,3,4,...
  #         BUILDS NESTED CFS VARIABLES BASED ON BEMT MATRIX
  #         REQUIRES OBJECT BEMT - THE CFS CLASSIFICATIONS
  #         Z1 IS THE EOSD CLASSIFICATION - USE ONLY WITH SINGLE
  #           Z VARIABLE AS IT IS NOT IN THE NESTED STRUCTURE
  #
  #############################################################
  
  cat("\n--- PREPARING QKEP MATRIX ---")
 
  cat("\nImage Categories: ", sort(unique(as.vector(BE))))
  
  # RESET THE FINAL MASK TO ENCOMPASS ALL
  FINALMASK <- (BE * 0) + 1
  
  # IF AUTOMASK IS ON - COMPUTE THE AUTOMASK
  if(AUTOMASK) {
    cat("\nUsing the auto-mask.")
    AMASK <- automask(IMG=BE)
  }
  else {
    cat("\nNot using the auto-mask.")
    AMASK <- (BE * 0) + 1
  }
  
  if(!is.null(GENERALMASK)) {
    # COMBINE THE AUTOMASK AND THE GENERAL MASK
    cat("\nCombining the general and auto-masks.")
    FINALMASK <- GENERALMASK * AMASK
  }
  else {
    cat("\nNo general mask.")
    FINALMASK <- AMASK
  }
  
  # AT THIS POINT, FINALMASK IS THE COMBINATION OF THE AUTO-MASK
  # AND GENERAL MASKS AND WILL BE THE ONE APPLIED
  
  
  # CONVERT IMAGES TO TABULAR FORM
  cat("\nConverting image to tabular form.")
  KI <- benya(BE)
  cat("\nConverting mask to tabular form.")
  MASK <- benci(FINALMASK)
  
  
  # --- X VARIABLE RECODE ---
  # IF REQUIRED, FILL X VARIABLES, STARTING WITH HIGHEST ONE FIRST
  # E.G., X3, X2, AND X1 (THE IMAGE NUMBER) LAST
  if(XRECODE) {
    cat("\nRecoding X Variables.")
    for(repx in NUMXRECODES:1) {
      KI <- cbind(XLUT[1,repx], KI)
    }
  }
  # FILL X1 COLUMN - IMAGE NUMBER
  cat("\nAdding X1 variable column - filled with 1s.")
  KI <- cbind(1,KI)
  
  
  # STORE TABULAR IMAGE DIMENSIONS
  NC <- dim(KI)[2]
  NR <- dim(KI)[1]
  
  
     
  # DRAW THE IMAGE AND MASK
  imaks(BE)
  #tkr(BE)
  title(sub="Original Image")
  #imaks(FINALMASK)
  imaks(FINALMASK)
  title(sub="Mask - Lt. Blue=Good")
  
  # APPLY THE MASK
  cat("\nApplying mask.")
  KI <- KI[MASK==1,]

    
  # STORE TABULAR IMAGE DIMENSIONS (AGAIN)
  NC <- dim(KI)[2]
  NR <- dim(KI)[1]
  
  
  # --- Z VARIABLE RECODE ---
  # IF REQUIRED, FILL Z VARIABLES, STARTING WITH LOWEST (COARSEST) ONE FIRST
  # E.G., Z2, Z3, ... (NOTE: Z1 IS ALWAYS THE ORIGINAL CLASSIFICATION FOR REFERENCE)
  if(ZRECODE) {
    cat("\nRecoding Z Variables.")
    # ADD THE NECESSARY NUMBER OF ZERO COLUMNS TO ACCOMODATE THE NEW Z VARIABLES
    for(newcol in 1:NUMZRECODES) {
      KI <- cbind(KI,0)
    }
  
    # ROW-BY-ROW, ADD PROPER LUT VALUES FOR NEW Z VARIABLES
    cat("\nApplying LUT for Z Variables.")
    for(row in 1:NR) {
      if(KI[row,NC] > 0) {
        KI[row, (NC+1):(NC+NUMZRECODES)] <- ZLUT[KI[row,NC],2:(NUMZRECODES+1)]
      }
    }
  }
  else {
    cat("\nNot recoding Z variables.")
    NUMZRECODES <- 0
  }

  # STORE TABULAR IMAGE DIMENSIONS
  NC <- dim(KI)[2]
  NR <- dim(KI)[1]
  
  
  # -----------------------------------------
  # IF WE WANT TO DROP THE ORIGINAL Z1 VARIABLE
  if(ZRECODE & DROP){
    # DROP Z1 TO REDUCE DIMENSIONALITY BUT SAVE NUMBER OF UNIQUE COLOURS FIRST
    assign(".COLOURS", max(KI[,NC-NUMZRECODES]), pos=1 )
    KI <- cbind(KI[,1:(NC - NUMZRECODES - 1)], KI[,(NC - NUMZRECODES + 1):NC])
    NUMZRECODES <- NUMZRECODES - 1
  }
  # -----------------------------------------
  
  
  # STORE TABULAR IMAGE DIMENSIONS
  NC <- dim(KI)[2]
  NR <- dim(KI)[1]
  
  # GENERATE NECESSARY VARIABLES FOR CREATING .QKEP
  assign(".N",dim(KI)[1],pos=1)
  # STORE THE NUMBER OF Y AND Z  VARIABLES
  if(ZRECODE) {
    if(XRECODE) {  # ZRECODE AND XRECODE
      if(DROP) {
        NNNX <- 1 + NUMXRECODES
        NNNY <- NC - NNNX - NUMZRECODES - 1
        NNNZ <- 1 + NUMZRECODES + 1
        STARTZ <- 2
      }
      else {
        NNNX <- 1 + NUMXRECODES
        NNNY <- NC - NNNX - (1 + NUMZRECODES)
        NNNZ <- 1 + NUMZRECODES
        STARTZ <- 1
      }
    }
    else {  # ZRECODE BUT NOT XRECODE
      if(DROP) {
        NNNX <- 1
        NNNY <- NC - NNNX - NUMZRECODES - 1 
        NNNZ <- 1 + NUMZRECODES + 1
        STARTZ <- 2
      }
      else {  # NO DROP
        NNNX <- 1
        NNNY <- NC - NNNX - (1 + NUMZRECODES)
        NNNZ <- 1 + NUMZRECODES
        STARTZ <- 1
      }
    }
  }
  else {  
    if(XRECODE) {  # XRECODE, NO ZRECODE
      NNNX <- 1 + NUMXRECODES
      NNNY <- NC - NNNX - 1
      NNNZ <- 1
      STARTZ <- 1
    }
    else {  # NO ZRECODE, NO XRECODE
      NNNX <- 1
      NNNY <- NC - NNNX - 1
      NNNZ <- 1
      STARTZ <- 1
    }
  }

  MAX <- apply(KI,2,max)
  KIKAP <- array(0,dim=MAX)
  KIKAP[KI] <- 1
  KIMOD <- dim(KIKAP)
  KIKAP <- KIKAP/sum(KIKAP)
  NAMESKI <- c(paste(sep="", "X", 1:NNNX), paste(sep="","Y",1:NNNY), paste(sep="", "Z", STARTZ:NNNZ) )
  cat(NAMESKI)
  names(KIMOD) <- NAMESKI

    # CREATE LOCAL VARIABLES
  assign(".QKEP", KIKAP, pos=1)
  assign(".LUT.Z", ZLUT, pos=1 )
  assign(".LUT.X", XLUT, pos=1 )
  if(!DROP) {
    assign(".COLOURS", max(KI[,NC-NUMZRECODES]), pos=1 )
  }
  assign(".KI", KI, pos=1)
  assign(".QND",KIMOD,pos=1)
  assign(".CIM",cim(BE),pos=1)
  return(cat("\n--- QKEP MATRIX BUILT ---\n"))
  
}

