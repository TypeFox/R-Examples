"bend4" <-
function (BE1=KEP1, BE2=KEP2, GENERALMASK=NULL, AUTOMASK=FALSE, XRECODE=FALSE, NUMXRECODES=1, XLUT=NULL, ZRECODE=FALSE, NUMZRECODES=1, ZLUT=NULL, DROP=TRUE, X1OMIT=FALSE) {

  #############################################################
  # 
  # TITLE:  t.bend4()
  # AUTHOR: TARMO REMMEL
  # DATE:   18 JAN, 2004
  # CALLS:  benya(), kateg(), benci(), t.automask()
  # NEEDS:  BEMT, INPUT IMAGES AS MATRIX OBJECTS, MASKS ARE OPTIONAL
  # NOTES:  READS TWO IMAGES INTO A MULTIDIMENSIONAL ARRAY
  #         COLOUR CODES MUST BE 1,2,3,4,...
  #         BUILDS NESTED CFS VARIABLES BASED ON BEMT MATRIX
  #         REQUIRES OBJECT BEMT - THE CFS CLASSIFICATIONS
  #         Z1 IS THE EOSD CLASSIFICATION - USE ONLY WITH SINGLE
  #           Z VARIABLE AS IT IS NOT IN THE NESTED STRUCTURE
  #
  #############################################################

  cat("\n--- PREPARING QKEP MATRIX ---")
  
  cat("\nImage 1 Categories: ", sort(unique(as.vector(BE1))))
  cat("\nImage 2 Categories: ", sort(unique(as.vector(BE2))))
  
  # RESET THE FINAL MASKS TO ENCOMPASS ALL
  FINALMASK1 <- (BE1 * 0) + 1
  FINALMASK2 <- (BE2 * 0) + 1
  
  # IF AUTOMASK IS ON - COMPUTE AUTO-MASK
  if(AUTOMASK) {
    cat("Using the auto-mask.")
    AMASK1 <- automask(IMG=BE1)
    AMASK2 <- automask(IMG=BE2)
  }
  else {
    cat("\nNot using the auto-mask")
    AMASK1 <- (BE1 * 0) + 1
    AMASK2 <- (BE2 * 0) + 1
  }
  
  #COMBINE THE AUTOMASK AND THE GENERAL MASK
  if(!is.null(GENERALMASK)) {
    cat("\nCombining the general and auto-masks.")
    FINALMASK1 <- GENERALMASK * AMASK1
    FINALMASK2 <- GENERALMASK * AMASK2
  }
  else {
    cat("\nNo general mask.")
    FINALMASK1 <- AMASK1
    FINALMASK2 <- AMASK2
  }
  
  # AT THIS POINT, FINALMASK IS THE COMBINATION OF THE AUTO-MASK
  # AND GENERAL MASKS AND WILL BE THE ONE APPLIED
  

  # CONVERT IMAGES TO TABULAR FORM
  cat("\nConverting images to tabular form.")  
  KI1 <- benya(BE1)
  KI2 <- benya(BE2)
  cat("\nConverting masks to tabular form.")  
  MASK1 <- benci(FINALMASK1)
  MASK2 <- benci(FINALMASK2)
  
  
  # --- X VARIABLE RECODE ---
  # IF REQUIRED, FILL X VARIABLES, STARTING WITH THE HIGHEST ONE FIRST
  # E.G., X3, X2, AND X1 (THE IMAGE NUMBER) LAST
  if(XRECODE) {
    cat("\nRecoding X Variables.")
    for(repx in NUMXRECODES:1) {
      KI1 <- cbind(XLUT[1,repx],KI1)
      KI2 <- cbind(XLUT[2,repx],KI2)
    }
  }
  # FILL X1 COLUMN - IMAGE NUMBER
  if(X1OMIT) {
    cat("\nFilling X1 variable column with 1s - OMITTING X1 TO REDUCE MEMORY USE")
    KI1 <- cbind(1,KI1)
    cat("\nFilling X1 variable column with 1s - OMITTING X1 TO REDUCE MEMORY USE")
    KI2 <- cbind(1,KI2)
  }
  else {
    cat("\nAdding X1 variable column - filled with image number.")  
    KI1 <- cbind(1, KI1)
    KI2 <- cbind(2, KI2)
  }  
  
  
  # STORE TABULAR IMAGE DIMENSIONS
  NC1 <- dim(KI1)[2]
  NR1 <- dim(KI1)[1]
  NC2 <- dim(KI2)[2]
  NR2 <- dim(KI2)[1]
  
  
  # DRAW THE IMAGES AND MASKS
  imaks(BE1)
  title(sub="Original Image 1")
  imaks(FINALMASK1)
  title(sub="Mask 1 - Lt. Blue=Good")
  imaks(BE2)
  title(sub="Original Image 2")
  imaks(FINALMASK2)
  title(sub="Mask 2 - Lt. Blue=Good")
  
  
  # APPLY THE MASKS
  cat("\nApplying masks.")
  KI1 <- KI1[MASK1==1,]
  KI2 <- KI2[MASK2==1,]
  
  
  # STORE TABULAR IMAGE DIMENSIONS AGAIN
  NC1 <- dim(KI1)[2]
  NR1 <- dim(KI1)[1]
  NC2 <- dim(KI2)[2]
  NR2 <- dim(KI2)[1]
  
  
  # --- Z VARIABLE RECODE ---
  # IF REQUIRED, FILL Z VARIABLES, STARTING WITH THE LOWEST (COARSEST) ONE FIRST
  # E.G., Z2, Z3, ... (NOTE: Z1 IS ALWAYS THE ORIGINAL CLASSIFICATION FOR REFERENCE)
  if(ZRECODE) {
    cat("\nRecoding Z Variables.")
    # ADD THE NECESSARY NUMBER OF ZERO COLUMNS TO ACCOMODATE THE NEW Z VARIABLES
    for(newcol in 1:NUMZRECODES) {
      KI1 <- cbind(KI1, 0)
      KI2 <- cbind(KI2, 0)
    }
    # ROW-BY-ROW, ADD PROPER LUT VALUES FOR NEW Z VARIABLES
    cat("\nApplying LUT for Z Variables.")
    for(row in 1:NR1) {
      if(KI1[row,NC1] > 0) {
        KI1[row, (NC1+1):(NC1+NUMZRECODES)] <- ZLUT[KI1[row,NC1], 2:(NUMZRECODES+1)]
      }
    }
    for(row in 1:NR2) {
      if(KI2[row,NC2] > 0) {
        KI2[row, (NC2+1):(NC2+NUMZRECODES)] <- ZLUT[KI2[row,NC2], 2:(NUMZRECODES+1)]
      }
    }   
  }
  else {
    cat("\nNot recoding Z Variables.")
    NUMZRECODES <- 0
  }

  
  # COMBINE IMAGE 1 AND 2 AND STORE DIMENSIONS
  KI <- rbind(KI1,KI2)
  NR <- dim(KI)[1]
  NC <- dim(KI)[2]

  
  # -----------------------------------------
  # IF WE WANT TO DROP THE ORIGINAL Z1 VARIABLE
  if(ZRECODE & DROP){
    # DROP Z1 TO REDUCE DIMENSIONALITY BUT SAVE UNIQUE NUMBER OF COLOURS FIRST
    assign(".COLOURS", max(KI[,NC-NUMZRECODES]), pos=1)
    KI <- cbind(KI[,1:(NC - NUMZRECODES - 1)], KI[,(NC - NUMZRECODES + 1):NC])
    NUMZRECODES <- NUMZRECODES - 1
  }
  # -----------------------------------------

  # STORE KI DIMENSIONS AGAIN
  NR <- dim(KI)[1]
  NC <- dim(KI)[2]
     
  
  # GENERATE NECESSARY VARIABLES FOR CREATING .QKEP
  assign(".N",dim(KI)[1],pos=1)
  # STORE THE NUMBER OF X, Y, AND Z VARIABLES
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
  
    
  names(KIMOD) <- NAMESKI
  
  # CREATE LOCAL VARIABLES
  assign(".QKEP",KIKAP, pos=1)
  assign(".LUT.Z", ZLUT, pos=1)
  assign(".LUT.X", XLUT, pos=1)
  if(!DROP) {
    assign(".COLOURS", max(KI[,NC-NUMZRECODES]), pos=1)
  }
  assign(".QND",KIMOD, pos=1)
  assign(".CIM",paste(sep="",cim(BE1),"\n",cim(BE2)),pos=1)
  return(cat("\n--- QKEP MATRIX BUILT ---\n"))

}

