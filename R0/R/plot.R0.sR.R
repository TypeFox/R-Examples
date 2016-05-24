# Name   : plot.sR
# Desc   : A tweaked "plot" function designed to easily plot all R objects from
#          a R0.sR-class list (result of est.R0).
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################

# Function declaration

plot.R0.sR <- function#Plot the R0/Rt value along with confidence interval of all requested models to epidemic data
### Plot the R0/Rt value along with confidence interval of all requested models to epidemic data
##details<< For internal use. Called by plot.
##keyword<< internal

(x, ##<<  Result of est.R0 (class sR)
 xscale="w", ##<< Scale to be adjusted on X axis. Can be "d" (day), "w" (week (default)), "f" (fornight), "m" (month).
 TD.split=FALSE, ##<< Parameter to force the display of both R(t) and the epidemic curve in the same window for TD method
 ... ##<< parameters passed to plot.R
) 
##details<< Tweaked plot() function that draws the reproduction number values for each method contained in the object constructed by est.RO().

  
# Code
  
{
	#Make sure x is of the right class.
	if (class(x)!="R0.sR") {
    stop("'x' must be of class 'R0.sR'")
	}
  
  #If invalid scale parameter is used, stop.
  if (xscale != "d" & xscale != "w" & xscale !="f" & xscale != "m") {
    stop("Invalid scale parameter.")
  }

  #Successive plots of individual model
  if (exists("EG", where = x$estimates)) {
    plot(x$estimates$EG, xscale=xscale, ...)
  }
  
  if (exists("ML", where = x$estimates)) {
    #x11()
    dev.new()
    plot(x$estimates$ML, xscale=xscale, ...)
  }
  
  if (exists("AR", where = x$estimates)) {
    #x11()
    dev.new()
    plot(x$estimates$AR, xscale=xscale, ...)
  }
  
  if (exists("TD", where = x$estimates)) {
    #x11()
    dev.new()
  plot(x$estimates$TD, xscale=xscale, TD.split=TD.split, ...)
  }
  
  if (exists("SB", where = x$estimates)) {
    #x11()
    dev.new()
  plot(x$estimates$SB, xscale=xscale, ...)
  }
  
  ### Called for its side effect :
  ### Draws all R0 or R(t) values from requested estimation methods.
}
