"bendManyX4" <-
function (INOBJ=XOBJ, GENERALMASK=NULL, AUTOMASK=F, XRECODE=F, NUMXRECODES=1, XLUT=NULL, ZRECODE=F, NUMZRECODES=1, ZLUT=NULL, DROP=T, X1OMIT=F) {
  
  #############################################################
  # 
  # TITLE:  t.bendManyX4()
  # AUTHOR: TARMO REMMEL
  # DATE:   28 OCT, 2004
  # CALLS:  importascii(), benya(), dups(), t.automask()
  # NEEDS:  
  # NOTES:  INOBJ = OBJ     AN OBJECT OF FILENAMES TO IMPORT AND PROCESS
  #			MASK			A MATRIX MASK OBJECT {0,1}
  #			USEMASK			BOOLEAN CONTROL FLAG TO USE MASK - IF MASK OBJECT IS NULL,
  #							  THEN, AN AUTOMASK IS GENERATED
  #         
  #         ONE IMAGE VARIABLE WITH MULTIPLE LEVELS (X=1,2,3,...)
  #         READS MULTIPLE IMAGES INTO A MULTIDIMENSIONAL ARRAY AND BUILDS .QKEP
  #         COLOUR CODES MUST BE 1,2,3,4,...
  #
  #############################################################
  
  cat("\n--- PREPARING QKEP MATRIX ---")
  
  # STORE THE NUMBER OF IMAGES IN INOBJ
  NX <- length(INOBJ)
  
  # LOOP FOR EACH IMAGE OBJECT
  for(rep in 1:NX) {

    cat("\nProcessing map ", rep, " of ", NX, sep="")
    
    # READ IMAGE OBJECT
    BE <- get(XOBJ[rep])
    
 
    # RESET THE FINAL MASK TO ENCOMPASS ALL
    FINALMASK <- (BE * 0) + 1
    
    # IF AUTOMASK IS ON - COMPUTE THE AUTOMASK
    if(AUTOMASK) {
      cat("\nUsing the auto-mask.")
      AMASK <- automask(IMG=BE,minvalue=1)
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
    # AND GENERAL MASK AND WILL BE THE ONE APPLIED
    
      
    # CONVERT IMAGES TO TABULAR FORM
    cat("\nConverting image to tabular form.")
    KItemp <- benya(BE)
    cat("\nConverting mask to tabular form.")
    MASK <- benci(FINALMASK)
 
   
    # --- X VARIABLE RECODE ---
    # IF REQUIRED, FILL X VARIABLES, STARTING WITH HIGHEST ONE FIRST
    # E.G., X3, X2, AND X1 (THE IMAGE NUMBER) LAST
    if(XRECODE) {
      cat("\nRecoding X Variables.")
      for(repx in NUMXRECODES:1) {
        KItemp <- cbind(XLUT[rep,repx], KItemp)
      }
    }

    if(X1OMIT) {
      cat("\nFilling X1 variable column with 1s - OMITTING X1 TO REDUCE MEMORY USE")
      KItemp <- cbind(1,KItemp)
    }
    else {
      # FILL X1 COLUMN - IMAGE NUMBER
      cat("\nAdding X1 variable column - filled with ", rep, "s.", sep="")
      KItemp <- cbind(rep,KItemp)
    }
  
    # STORE TABULAR IMAGE DIMENSIONS
    NC <- dim(KItemp)[2]
    NR <- dim(KItemp)[1]
  
    
    # DRAW THE IMAGE AND MASK
    imaks(BE)
    title(sub="Original Image")
    imaks(FINALMASK)
    title(sub="Mask - Lt. Blue=Good")
  
    # APPLY THE MASK
    cat("\nApplying mask.")
    KItemp <- KItemp[MASK == 1,]

    # STORE TABULAR IMAGE DIMENSIONS (AGAIN)
    NC <- dim(KItemp)[2]
    NR <- dim(KItemp)[1]
  

    cat("\nImage Ctegories:", sort(unique(as.vector(KItemp[,NC]))) )

    # --- Z VARIABLE RECODE ---
    # IF REQUIRED, FILL Z VARIABLES, STARTING WITH LOWEST (COARSEST) ONE FIRST
    # E.G., Z2, Z3, ... (NOTE: Z1 IS ALWAYS THE ORIGINAL CLASSIFICATION FOR REFERENCE)
    if(ZRECODE) {
      cat("\nRecoding Z Variables.")
      # ADD THE NECESSARY NUMBER OF ZERO COLUMNS TO ACCOMODATE THE NEW Z VARIABLES
      for(newcol in 1:NUMZRECODES) {
        KItemp <- cbind(KItemp,0)
      }
      # ROW-BY-ROW, ADD PROPER LUT VALUES FOR NEW Z VARIABLES
      cat("\nApplying LUT for Z Variables.")
      for(row in 1:NR) {
        if(KItemp[row,NC] > 0) {
          KItemp[row, (NC+1):(NC+NUMZRECODES)] <- ZLUT[KItemp[row,NC],2:(NUMZRECODES+1)]
        }
      }
    }
    else {
      cat("\nNot recoding Z variables.")
      NUMZRECODES <- 0
    }
  
          
    # ADD IMAGE TO MASTER MATRIX
    if(rep == 1) {
      # IF FIRST IMAGE
      cat("\nAdding first image to the master matrix.")
      KI <- KItemp
    }
    else {
      # IF SUBSEQUENT IMAGES
      cat("\nAdding image ",rep, " to the master matrix.", sep="")
      KI <- rbind(KI,KItemp)
    }

    cat("\n")
    
  } # END FOR
  
  
  # STORE MASTER MATRIX DIMENSIONS
  NR <- dim(KI)[1]       
  NC <- dim(KI)[2]
  
  
  # -----------------------------------------
  # IF WE WANT TO DROP THE ORIGINAL Z1 VARIABLE
  if(ZRECODE & DROP){
    # DROP Z1 TO REDUCE DIMENSIONALITY BUT SAVE THE NUMBER OF UNIQUE COLOURS FIRST
    assign(".COLOURS", max(KI[,NC-NUMZRECODES]), pos=1)
    KI <- cbind(KI[,1:(NC - NUMZRECODES - 1)], KI[,(NC - NUMZRECODES + 1):NC])
    NUMZRECODES <- NUMZRECODES - 1
  }
  # -----------------------------------------

  # STORE MASTER MATRIX DIMENSIONS
  NR <- dim(KI)[1]       
  NC <- dim(KI)[2]
  
  
    
  # GENERATE NECESARY VARIABLES FOR CREATING .QKEP
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

  if(X1OMIT) {
    # COMPENSATE FOR DUPLICATE ROW ENTRIES IF X1 IS OMITTED
    dups(KI)
    MAX <- apply(.KIUNIQUE, 2, max)
    KIKAP <- array(0, dim=MAX)
    KIKAP[.KIUNIQUE] <- .COUNT
  }
  else {
    # IF X1 IS NOT OMITTED
    MAX <- apply(KI,2,max)
    KIKAP <- array(0,dim=MAX)
    KIKAP[KI] <- 1 
  }
  
  KIMOD <- dim(KIKAP)
  KIKAP <- KIKAP/sum(KIKAP)
  NAMESKI <- c(paste("X",1:NNNX,sep=""), paste("Y", 1:NNNY, sep=""), paste("Z", STARTZ:NNNZ, sep="") )
  names(KIMOD) <- NAMESKI

  
  # BUILD .CIM OBJECT CONTAINING OBJECT NAMES FOR EACH INPUT IMAGE
  tempCIM <- NULL
  for(i in 1:NX) {
    tempCIM <- paste(tempCIM, "\n", cim(get(XOBJ[i])), sep="")
  }
  tempCIM <- substring(tempCIM, 2, nchar(tempCIM))
  
      
  # MAKE LOCAL COPIES FOR LATER USE
  assign(".QKEP",KIKAP, pos=1)
  assign(".LUT.Z", ZLUT, pos=1)
  assign(".LUT.X", XLUT, pos=1)
  if(!DROP) {
    assign(".COLOURS", max(KI[,NC-NUMZRECODES]), pos=1)
  }
  assign(".QND",unlist(KIMOD), pos=1)
  #assign(".CIM",paste(sep="",cim(XOBJ), pos=1))
  assign(".CIM", tempCIM, pos=1)
  return(cat("\n--- QKEP MATRIX BUILT ---\n"))
}

