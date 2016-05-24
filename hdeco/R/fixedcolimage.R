"fixedcolimage" <-
function (BE=checker(8)) { 

  #############################################################
  # 
  # TITLE:  fixedcolimage()
  # AUTHOR: TARMO REMMEL (MODIFIED FROM: SANDOR KABOS) 
  # DATE:   28 OCT 2004
  # CALLS:  N/A
  # NEEDS:  IMAGE MATRIX OBJECT
  # NOTES:  PROVIDES PROPER ORIENTATION OF IMAGE FOR VIEWING
  #         USES A CONSISTENT LEGEND
  #         BASED ON imaks()
  #
  #############################################################
  
  LENG <- 4
  szinszam <- 6
  
  "sorrev" <- function (BE=.Q) { 
    NR<-dim(BE)[1] 
    NRR<-1:NR 
    REVENR<-rev(NRR) 
    BE[REVENR,] 
  } 

  NR <- dim(BE)[1] 
  NC <- dim(BE)[2] 
  LENG <- min(LENG,NR,NC) 
  par(pty="s") 
  LR <- NR%/%LENG 
  LNR <- 1:LENG 
  LNR <- LNR * LR 
  LNR <- c(1, LNR) 
  LABR <- 1:NR 
  LABR <- LABR[LNR] 
  LABR <- as.character(LABR) 
  ATR <- seq(0, 1, length=(NR)) 
  ATR <- ATR[LNR] 
  
  # OBTAIN THE UNIQUE COLORS FROM THE IMAGE (MAP)
  uniquecolors <- as.vector(as.integer(levels(factor(BE))))
  cat("\nOriginal Image Categories: ", uniquecolors, sep=", ")

  # CHECK FOR ONLY ONE COLOUR (PALETTES NEED A MINIMUM OF 2 COLOURS)
  if(length(uniquecolors) == 1) {
    cat("\nOnly one colour in this map.\n")
  }
    
  # BUMP COLOURS UP BY ONE SUCH THAT 0 (ZERO) IS 1 IF NECESSARY, AND WILL BE DRAWN IN WHITE (BACKGROUND)
  if(uniquecolors[1] == 0) {
    cat("\nNoData Values exist in map.")
    uniquecolors <- uniquecolors + 1
    BE <- BE + 1
    # SET THE CUSTOM PALETTE COLOURS FOR CATEGORIES 1, 2, ... 
    custompalette <- colors()[c(1,490,32,254,81,652,547,646,610,625,536,508,530,498,479,410,181,136,142,20,11)]
    cat("\nStarting Palette: ", custompalette, sep=", ")
  }
  else {
    # SET THE CUSTOM PALETTE COLOURS FOR CATEGORIES 1, 2, ... 
    custompalette <- colors()[c(490,32,254,81,652,547,646,610,625,536,508,530,498,479,410,181,136,142,20,11)]
    cat("\nStarting Palette: ", custompalette, sep=", ")
  }
  
  # SUBSET THE CUSTOM PALETTE FOR ONLY THOSE COLOURS USED IN THE MAP BEING DISPLAYED]
  custompalette <- custompalette[uniquecolors]
  cat("\nExtracted Palette: ", custompalette, sep=", ")    
   
  # SET THE SUBSET PALETTE AS THE CURRENT PALETTE
  palette(custompalette)

  imgcolours <- as.integer(levels(factor(BE)))
  for(a in 0:length(imgcolours)){
    BE[BE == imgcolours[a]] <- a
  }
    
  # DRAW THE IMAGE
  image(z=t(sorrev(BE)), col=palette(), xaxt="n", yaxt="n", axes=F) 
  axis(1, at=ATR, labels=LABR) 
  axis(2, at=(1-ATR), labels=(LABR)) 
  title(cim(BE)) 
  par(pty="m")
  cat("\nDone.")
   
}

