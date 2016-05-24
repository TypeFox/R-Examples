"checker" <-
function(block=4) {

  #############################################################
  # 
  # TITLE:  checker()
  # AUTHOR: TARMO REMMEL
  # DATE:   11 FEB 2004
  # CALLS:  fixedcolimage()
  # NEEDS:  
  # NOTES:  CREATES CHECKERBOARD PATTERN IMAGES BASED ON BLOCK
  #		SIZE - ALL REUSLTS ARE BINARY
  #
  #############################################################

  # BUILDS 64x64 CHECKER BOARDS

  IMG <- NULL
  
  par(mfrow=c(2,2))
    
  black <- rep(1:2,(64/block/2), each=block)
  dim(black) <- c(64,1)
  white <- rep(2:1, (64/block/2), each=block)
  dim(white) <- c(64,1)
  
  colour <- black
  for(rmain in 1:(64/block)) {
    for(r in 1:block) {
      IMG <- cbind(IMG,colour)
    }
    
    if(colour == black) {
      colour <- white
    }
    else {
      colour <- black
    }
   
  }
  
  fixedcolimage(IMG)
  return(IMG) 
 
}

