##  PLOT BACKGROUND FUNCTIONS
##  Kyle Elmy and Jim Kaklamanos
##  1 December 2015

##  NOTE:
##  This file contains functions that are used to create the logarithmic
##  and semilogarithmic plots that are used by other functions in the package.


########################################################################################
##  1. BACKGROUND FUNCTION logAxis.plot = plotting function (one axis at a time)
##
##  Arguments:  data = vector of data values for plot;
##              type = "x" for horizontal axis; "y" for vertical axis
##              vertical = TRUE for vertical gridlines across plot;
##                         FALSE otherwise (default = FALSE)

logAxis.plot <- function(data, type, gridlines){
  
  ##  Range of data
  pow.min <- 10^floor(log10(min(data[data > 0])))
  pow.max <- 10^ceiling(log10(max(data)))
  
  ##  Tick marks
  major.ticks <- 10^seq(from = log10(pow.min), to = log10(pow.max), by = 1)
  minor.ticks <- c()
  j <- 1
  for(i in 1:(length(major.ticks)-1)){
    minor.ticks[j:(j+9)] <- seq(from = major.ticks[i], to = major.ticks[i+1],
                                by = major.ticks[i])
    j <- j + 10
  }
  labels <- formatC(major.ticks, format = "fg")

  ##  Side for axis
  if(type == "x"){
    side <- 1
  } else{
    if(type == "y"){
      side <- 2
    }
  }
  
  ##  Draw axes
  axis(side = side, at = major.ticks, labels = labels, tcl = -0.6)
  axis(side = side, at = minor.ticks, labels = FALSE, tcl = -0.3)

  ##  Draw gridlines
  if(gridlines == TRUE){
    axis(side = side, at = major.ticks, labels = FALSE, tck = 1,
         col = "gray70", lwd = 2)
    axis(side = side, at = minor.ticks, labels = FALSE, tck = 1,
         col = "gray50", lwd = 1)
  }

}




########################################################################################
##  2. BACKGROUND FUNCTION logAxis = wrapper function for performing calculations
##
##  Arguments:  x = data for x-axis, y = data for y-axis
##              gridX and gridY are TRUE if gridlines in the vertical or horizontal
##                  directions are to be drawn; FALSE otherwise
##              

logAxis <- function(x = NA, y = NA, gridX = FALSE, gridY = FALSE){

  if(is.na(x) == FALSE && is.na(y) == TRUE){
    logAxis.plot(data = x, type = "x", gridlines = gridY)
    
  } else{
    if(is.na(x) == TRUE && is.na(y) == FALSE){
      logAxis.plot(data = y, type = "y", gridlines = gridX)
      
    } else{
      if(is.na(x) == FALSE && is.na(y) == FALSE){
        logAxis.plot(data = x, type = "x", gridlines = gridY)
        logAxis.plot(data = y, type = "y", gridlines = gridX)
      }
    }
  }
}
