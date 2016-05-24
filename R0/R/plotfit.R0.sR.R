# Name   : plotfit.R0.sR
# Desc   : A tweaked "plot" function designed to easily plot all R objects from
#          a sR-class list (result of est.R0).
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################

# Function declaration

plotfit.R0.sR <- function#Plot the fit of all requested models to epidemic data
### Plots the fit of all requested models to epidemic data
##details<< For internal use. Called by plotfit.
##keyword<< internal

(x, ##<< Result of est.R (class R0.R)
 all=TRUE, ##<< Should the whole epidemic curve be shown
 xscale="w", ##<< Scale to be adjusted on X axis. Can be "d" (day), "w" (week (default)), "f" (fornight), "m" (month).
 SB.dist=TRUE, ##<< Should R distribution throughout the epidemic be plotted for SB method? (default: TRUE)
 ... ##<< Parameters passed to plot
 ) 
  ##details<< Tweaked plot() function that draws the epidemic data and model fit of each method contained in the object constructed by est.RO().
  
  
  # Code
  
{
  #Make sure x is of the right class.
  if (class(x)!="R0.sR") {
    stop("'x' must be of class 'sR'")
  }
  
  #If invalid scale parameter is used, stop.
  if (xscale != "d" & xscale != "w" & xscale !="f" & xscale != "m") {
    stop("Invalid scale parameter.")
  }
  
  #Successive plots of individual model
  if (exists("EG", where = x$estimates)) {
    plotfit(x$estimates$EG, xscale=xscale, ...)
  }
  
  if (exists("ML", where = x$estimates)) {
    #x11()
    dev.new()
    plotfit(x$estimates$ML, xscale=xscale, ...)
  }
  
  if (exists("TD", where = x$estimates)) {
    #x11()
    dev.new()
    plotfit(x$estimates$TD, xscale=xscale, ...)
  }
  
  if (exists("SB", where = x$estimates)) {
    #x11()
    dev.new()
    plotfit(x$estimates$SB, xscale=xscale, SB.dist = SB.dist, ...)
  }
  
  ### Called for its side effect :
  ### Draws all the epidemic curves and associated fit data computed by est.R0
}
