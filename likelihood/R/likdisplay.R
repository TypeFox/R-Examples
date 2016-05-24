#######################################################
# likdisplay
# Controls the output display of the annealing run. This creates new plots
# each time it is called.  It updates both the graph showing the likelihood
# curve and the display with the current run stats.
#
# Author:  Lora Murphy, Cary Institute of Ecosystem Studies
# murphyl@caryinstitute.org
#######################################################
likdisplay<-function(lhistcycles, lhisthood, slp, R2, aiccorr, temp, max_cycles) {

  # If the user has closed the window, don't graph
  if(length(.Devices) == 1 || dev.cur() == 1) {
    return()
  }

  # Plot the likelihood history
  ylim<-c(max(min(lhisthood), max(lhisthood)*2), max(lhisthood)/1.1)
  # Make the y axis such that the minimum is no less than 50% of the maximum
  plot(lhistcycles, lhisthood, type="o", main="Likelihood History",
       xlab="Annealing Cycles", ylab="Log Likelihood", 
       ylim=c(max(min(lhisthood), max(lhisthood)*2), max(lhisthood)/1.1))

  # Plot the text - figure out the ymax
  # Find out how many lines of text we can fit in the display
  ycoord <- par("csi")
  numlines <- floor(par("pin")[[2]] / ycoord)
  ymax<-numlines

  # Create an empty plot
  plot(0,0,type="n", axes=FALSE, xlab="",ylab="", 
       xlim=c(0,5), ylim=c(0,ymax))
      
  # Write out the text - only as many lines as will fit
  if (5 < numlines)
    text(0.1,5, paste("Completed cycle ", formatC(lhistcycles[[length(lhistcycles)]], format="f", digits=2), " of ", formatC(max_cycles, format="f", digits=0), ""), pos=4)
  if (4 < numlines)
    text(0.1,4,paste("Current temp: ", formatC(temp, format="g", digits=2), ""), pos=4)
  if (3 < numlines)
    text(0.1,3,paste("Current best likelihood: ", round(lhisthood[[length(lhisthood)]], 2), ""), pos=4)
  if (2 < numlines)
    text(0.1,2,"Goodness of fit - regression of observed on predicted:", pos=4)
  if (1 < numlines)
    text(0.1,1,paste("Slope:  ", formatC(slp, format="g", digits=2), "  R2:  ", formatC(R2, format="g", digits=2), ""), pos=4)
  if (0 < numlines)
    #This next one's Y position is not 0 because sometimes it cut off the bottom of the letters
    text(0.1,0.1,paste("AIC corr: ", formatC(aiccorr, format="f", digits=2), ""), pos=4)
}