# Name   : plotfit.R
# Desc   : S3 Method for plotting the fit of an epidemic model
# Date   : 2012/05/02
# Author : Boelle, Obadia
###############################################################################


# Function declaration

plotfit <- function#Generic S3 method to plot either "R0.R" and "R0.sR" objects
###Generic S3 method to plot either "R0.R" and "R0.sR" objects

(x, ##<< Object for which the fit should be plotted.
 all=TRUE, ##<< Should the whole epidemic curve be shown
 xscale="w", ##<< Scale to be adjusted on X axis. Can be "d" (day), "w" (week (default)), "f" (fornight), "m" (month).
 SB.dist=TRUE, ##<< Should R distribution throughout the epidemic be plotted for SB method? (default: TRUE)
 ... ##<< parameters passed to plot.R
 ) 
  ##details<< plot.fit is designed to either call plot.fit.R0.R or plot.fit.R0.sR.
  ## This S3 Method allows for plottinf the goodness of fit of a model to the original epidemic curve provided by user.
  ## Depending on the method of estimation, the graphical output will vary:
  ## - EG, ML and TD methods will show the original epidemic curve, along with the best-fitting prediction model
  ## - AR will only show the epidemic curve, since no actual model is computed
  ## - RTB will display 9 density curves for the R distribution throughout the epidemic
  
  
  # Code
  
{
  UseMethod("plotfit")
}
