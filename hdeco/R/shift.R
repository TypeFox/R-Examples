"shift" <-
function(map=demoimage1, dir=1, n=1, draw=TRUE, verbose=TRUE) {

  #########################################################################
  #
  # TITLE:     shift()
  # AUTHOR:    TARMO REMMEL
  # DATE:      26 APRIL 2007
  # CALLS:     fixedcolimage()
  # NEEDS:     A CATEGORICAL MAP WITH EQUAL X AND Y DIMENSIONS AS A MATRIX
  # NOTES:     THE PARAMETER n DOMINATES dir; THUS IF n=0 BUT dir IS OUT
  #            OF RANGE, THE SCRIPT WILL STILL PRODUCE A VALID UNSHIFTED
  #            MAP
  #
  #            SET THE dir PARAMETER ACCORDING TO THE FOLLOWING:
  #            1 = UP (NORTH) 
  #            2 = RIGHT (EAST)
  #            3 = DOWN (SOUTH)
  #            4 = LEFT (WEST)
  #
  # SAMPLE:    shift(demoimage1, dir=2, n=3, draw=TRUE, verbose=TRUE)
  #
  #########################################################################
  
  # INITIALIZE INFORMATION VECTOR FOR TEXTUAL OUTPUT
  direction <- c("up", "right", "down", "left")


  # HANDLE n=0 (NO SHIFT SPECIFIED)
  if(n == 0) {
    if(verbose) {
      cat("\nNo shift specified.  Output is the same as input.\n")
    }
    shiftmap <- map
    par(mfrow=c(1,2))
    if(draw) { fixedcolimage(map) }
    attr(shiftmap, "cim") <- paste("Shifted Map: ", n, " ", direction[dir], sep="")
    if(draw) { fixedcolimage(shiftmap) }
    return(shiftmap)
  }
 
  # PROCESS ONLY IF CORRECT DIRECTIONAL INTEGER IS PROVIDED
  if( (dir %in% c(1,2,3,4)) & (n >= 0) & (n < dim(map)[1]) ) {

    if(verbose) {
      cat("\nWrapping map ", n, " unit(s) ", direction[dir], ".\n", sep="")
    }
  
    # PREPARE DRAWING ENVIRONMENT AND DRAW THE ORIGINAL IMAGE
    par(mfrow=c(1,2))
    if(draw) { fixedcolimage(map) }

    # SHIFT UP
    if(dir == 1) {
      shiftmap <- rbind(map[(n+1):dim(map)[1],], map[1:n,])
    }

    # SHIFT RIGHT
    if(dir == 2) {
      shiftmap <- cbind(map[,((dim(map)[2]+1)-n):(dim(map)[2])], map[,1:(dim(map)[2]-n)])
    }

    # SHIFT DOWN
    if(dir == 3) {
      shiftmap <- rbind(map[((dim(map)[1]+1)-n):(dim(map)[1]),], map[1:(dim(map)[1]-n),])
    }

    # SHIFT LEFT
    if(dir == 4) {
      shiftmap <- cbind(map[,(n+1):dim(map)[2]], map[,1:n])
    }

    # DRAW SHIFTED MAP WITH MEANINGFUL TITLE
    attr(shiftmap, "cim") <- paste("Shifted Map: ", n, " ", direction[dir], sep="")
    if(draw) { fixedcolimage(shiftmap) }
    return(shiftmap)
  }
  else {
    # PROVIDE FEEDBACK REGARDING BAD PARAMETER SETTINGS
    cat("\n PARAMETER ENTRY ERROR:")
    cat("\n Ensure that the dir parameter is either 1, 2, 3, or 4")
    cat("\n Ensure that the n parameter is > 0")
    cat("\n Ensure that the n parameter is < the dimensions of the input map\n\n")
    return(NULL)

  } # END IF

} # END FUNCTION

