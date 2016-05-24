#' Phase plot of BCPA output
#' 
#' Behavioral phase plot of BCPA output.  Mean and standard deviation are on the x and y axis.  Size and color of points represent the time scale of autocorrelation

#' @param windowsweep \code{windowsweep} object, i.e. the output of the \code{\link{WindowSweep}} function.
#' @param type whether to plot smooth or flat bcpa output
#' @param clusterwidth for flat BCPA, this is the temporal range within which change points are considered to be within the same cluster. 
#' @param ... additional arguments passed to the \code{plot} base function.
#' @param legend.where where to place the tau legend (see \code{\link{legend}}).
#' @seealso \code{\link{WindowSweep}}, \code{\link{PartitionParameters}}
#' @author Eliezer Gurarie
#' @examples
#' if(!exists("Simp.ws"))
#' {
#'  data(Simp)
#'  Simp.ws <- WindowSweep(GetVT(Simp), "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
#' }
#' 
#' PhasePlot(Simp.ws)

PhasePlot <- function(windowsweep, type = c("smooth","flat")[1], clusterwidth = 1, ..., legend.where="bottomright")
{
  partition <- PartitionParameters(windowsweep, type, clusterwidth)
  mu.hat <- partition$mu.hat
  s.hat <- partition$s.hat
  rho.hat <- partition$rho.hat
  
  if(type=="smooth")
  {
    rho.int <- round((rho.hat - min(rho.hat, na.rm=1))/(max(rho.hat, na.rm=1) - min(rho.hat, na.rm=1))*999+1)
    palette(topo.colors(1000, alpha=0.5) )
    
    plot(mu.hat, s.hat, type="l", col="grey", xlab=expression(hat(mu)),  ylab=expression(hat(sigma)), ...)
    points(mu.hat, s.hat, cex=sqrt((rho.hat - min(rho.hat, na.rm=1))/(max(rho.hat, na.rm=1) - min(rho.hat, na.rm=1))*4), 
           col=rho.int, pch=19)
    legend(legend.where, pt.cex = sqrt(c(.5,1,2,4)), pch = 19, title=expression(hat(tau)), col = 250*c(.5,1,2,4), bty="n",
           legend=signif(seq(.2,1,length=4)*max(rho.hat, na.rm=TRUE),2))
    palette("default")
    points(mu.hat[1], s.hat[1], pch=24, bg="darkgreen", cex=2, col="darkgrey")
    points(mu.hat[length(mu.hat)], s.hat[length(mu.hat)], pch=23, bg="red", cex=2,  col="darkgrey")
  }
  
  
  if(type=="flat")
  {
    rho.int <- round((rho.hat - min(rho.hat, na.rm=1))/(max(rho.hat, na.rm=1) - min(rho.hat, na.rm=1))*999+1)
    palette(topo.colors(1000, alpha=0.5) )
    
    N.reps <- as.vector(table(paste(mu.hat, s.hat, rho.hat)))
    N.count <- rep(N.reps, N.reps)
    
    plot(mu.hat, s.hat, type="l", col="grey", xlab=expression(hat(mu)),  ylab=expression(hat(sigma)), ...)
    points(mu.hat, s.hat, cex=sqrt(N.count), 
           col=rho.int, pch=19)
    legend(legend.where, pt.cex =2, pch = 19, title=expression(hat(tau)), col = 250*c(.5,1,2,4), 
           legend=signif(seq(.2,1,length=4)*max(rho.hat, na.rm=TRUE),2))
    palette("default")
  }
}