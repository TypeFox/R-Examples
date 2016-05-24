"imaks" <-
function (BE=demoimage1, ncolours=NULL, LENG=4) { 

  #############################################################
  # 
  # TITLE:  imaks()
  # AUTHOR: SANDOR KABOS, MODIFIED BY TARMO REMMEL
  # DATE:   15 AUG, 2003
  # CALLS:  N/A
  # NEEDS:  IMAGE MATRIX OBJECT
  # NOTES:  PROVIDES PROPER ORIENTATION OF IMAGE FOR VIEWING
  #         LENG CONTROLS THE LABELING
  #
  #############################################################
  
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
  
  # IF SZINSZAM IS NOT GIVEN IN FUNCTION CALL, SET IT TO THE MAXIMUM NUMBER OF COLOURS IN THE IMAGE
  if(is.null(ncolours)) {
    ncolours <- length(table(BE))
  }
  
  # SET THE DEFAULT COLOUR SCHEME
  COOL <- topo.colors(ncolours) 
  
  # SET SPECIFIC COLOUR SCHEMES FOR 2-6 COLOUR IMAGES
  if(ncolours == 1) {
     COOL <- topo.colors(12)[c(4)]  # ENSURES FULL GOOD MASK IS Lt. BLUE - USUALLY NOT USED
  }
  else {
    if(ncolours == 2) {
      COOL <- topo.colors(12)[c(2,4)]
    }
    else {
      if(ncolours == 3) {
        COOL<- topo.colors(12)[c(2,4,6)]
      }
      else {
        if(ncolours == 4) {
          COOL <- topo.colors(12)[c(2,4,6,8)] 
        }
        else {
          if(ncolours == 5) {
            COOL <- topo.colors(12)[c(2,4,6,8,10)] 
          }
          else {
            if(ncolours == 6) {
              COOL <- topo.colors(12)[c(2,4,6,8,10,12)] 
            }
          }
        }
      }
    }
  }

     
  image(t(sorrev(BE)), col=COOL, xaxt="n", yaxt="n", axes=F) 
  axis(1, at=ATR, labels=LABR) 
  axis(2, at=(1-ATR), labels=(LABR)) 
  title(cim(BE)) 
  par(pty="m") 
  
}

