################################################################################
# Automation and Robotic Section
# European Space Research and Technology Center (ESTEC)
# European Space Agency (ESA)
##
# Centre for Automation and Robotics (CAR)
# CSIC- Universidad Politecnica de Madrid
#
# Header: plotAV(param)
#
# Author: Javier Hidalgo Carrio
#
# Date: 10-05-2010
#
# Input Parameters:
#       one data frame structure with tree fields
#         $av -> Allan Variance
#         $time -> cluster time of the computation
#         $error -> error of the estimation (quality of the variance)
#
# Output Parameters:
#       visual plotting
#
# Description:
# The R packages required for the execution of this script are:
#  gplots, install.packages("gplots", dependencies = TRUE), for installing
#
# plotAV function computes the plot of the values for a Allan variance 
# estimation.
#
# License: GPL-2
#
#' @export
#' @importFrom graphics plot lines axis grid title
#' @importFrom gplots plotCI
#
################################################################################
plotav <- function (avdf)
{
  #Load library
  #library(gplots)

  #Plot Error Bars and Confidence Intervals
  plotCI(x= avdf$time, y=sqrt(avdf$av), uiw =  sqrt(avdf$av)*avdf$error, gap=0, lty = 1, xaxt="n", yaxt="n", col="blue", log="xy", lwd=1,  xlab="", ylab="")
  #Plot the line
  lines (avdf$time,sqrt(avdf$av), col="blue")
  #Axis numbers scale
  axis(1, c(0.001,0.01, 0.1, 0, 1, 10, 100, 1000, 10000))
  axis(2, c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0, 1, 10, 100, 1000))
  #Grid log-log
  grid(equilogs=TRUE, lwd=1, col="orange")
  #Information (title and axis labels)
  title(main = "Allan variance Analysis", xlab = "Cluster Times (Sec)", ylab = "Allan Standard Deviation")

  #legend(10, 5e-03, c("Allan Variance"),  fill = c("blue"))
}
