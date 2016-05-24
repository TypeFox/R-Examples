#' Path plot of BCPA output
#' 
#' Plots the animal's trajectory, with segments color-coded according to the time scale / auto-correlation of the BCPA output, and width of segments proportional to the estimated mean of the BCPA.  
#' 
#' @param Data the track data to be plotted - most typically, output of the \code{\link{GetVT}} function.
#' @param windowsweep \code{windowsweep} object, i.e. the output of the \code{\link{WindowSweep}} function.
#' @param type whether to plot smooth or flat bcpa output
#' @param clusterwidth for flat BCPA, this is the temporal range within which change points are considered to be within the same cluster. 
#' @param plotlegend whether to plot a legend.
#' @param tauwhere where to place the legend for the time-scale / auto-correlation.  Can be one of "nowhere", "top", "bottom", "left", "right", "topleft", "topright", "bottomright", "bottomleft".
#' @param n.legend number of labels in legend. 
#' @param ncol.legend number of columns in the legend.
#' @param bty.legend whether to draw a box around the legend.
#' @param ... additional arguments to pass to the \code{plot} base function.
#' @author Eliezer Gurarie
#' @examples
#' if(!exists("Simp.ws"))
#' {
#'  data(Simp)
#'  Simp.ws <- WindowSweep(GetVT(Simp), "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
#' }
#' 
#' PathPlot(Simp, Simp.ws, plotlegend=TRUE, n.legend=3)
#' PathPlot(Simp, Simp.ws, type="flat", clusterwidth=3, plotlegend=TRUE)

PathPlot <-  function (Data, windowsweep, type=c("smooth", "flat")[1], clusterwidth = 1, plotlegend = FALSE, tauwhere = "topleft", n.legend = 5, ncol.legend = 1, bty.legend="n", ...) 
{
  if(!("Z" %in% names(Data)))
    z <- Data$X + 1i*Data$Y else z <- Data$Z
  
  if(type=="flat")
    pp <- PartitionParameters(windowsweep, type = type, clusterwidth=clusterwidth)
  
  if(type == "smooth")  
    if("pp.smooth" %in% names(windowsweep)) pp <- windowsweep$pp.smooth else pp <- PartitionParameters(windowsweep, type = type)
  
  Segments <- function(z, col=col, lwd=lwd)
  {
    n <- length(z)
    segments(Re(z[-n]), Im(z[-n]), Re(z[-1]), Im(z[-1]), col=col, lwd=lwd)
  }  
  
  mu.hat <- pp$mu.hat
  rho.hat <- pp$rho.hat
  
  rho.max <- max(rho.hat, na.rm=1)
  rho.int <- round(rho.hat/rho.max * 999 + 1)
  
  palette(topo.colors(1000))  
  plot(z, asp=1, pch=19, cex=0.5, col="grey", ...)
  points(z[c(1,length(z))], pch=c(24, 23), bg=c("green", "red"), cex=2, lwd=1.5, col="darkgrey")
  
  Segments(z, col=rho.int, lwd=abs(mu.hat/max(mu.hat, na.rm=TRUE))*4)
  
  if(plotlegend)
    legend(tauwhere,  lty=1, 
           title=expression(hat(tau)), ncol = ncol.legend, bty=bty.legend, lwd=2, 
           col=seq(0,999, length=n.legend) + 1, 
           legend = signif(seq(0, max(rho.hat), length = n.legend),2))
  palette("default")
}
