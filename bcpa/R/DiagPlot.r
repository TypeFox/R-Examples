#' Diagnostic plot for BCPA
#' 
#' Draws diagnostic plots for BCPA analysis.  Specifically: a qqplot, a histogram (with a N(0,1) density curve), and an acf of the standardized residuals of an analysis. 
#' 
#' @param windowsweep a \code{windowsweep} object, i.e. the output of the \code{\link{WindowSweep}} analysis.
#' @param type whether to diagnose the model fitted for a smooth or flat BCPA. 
#' @param plotme logical - whether or not to plot the diagnostics
#' @param values logical - whether or not to return the values of the standardized residuals.
#' @param ... additional arguments to pass to the \code{\link{PartitionParameters}} function.
#' @return If \code{values} is TRUE, returnn the values of the standardized residuals.
#' @seealso \code{\link{PartitionParameters}}
#' @author Eliezer Gurarie
#' @examples
#' data(Simp)
#' if(!exists("Simp.VT"))
#'  Simp.VT <- GetVT(Simp)
#' if(!exists("Simp.ws"))
#'  Simp.ws <- WindowSweep(Simp.VT, "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
#' DiagPlot(Simp.ws)
#' DiagPlot(Simp.ws, type="flat")
#' # The Simp's diagnostic plots are excellent.


DiagPlot <- function(windowsweep, type = c("smooth", "flat")[1], plotme=TRUE, values = FALSE, ...)
{
  x <- windowsweep$x
  pp <- PartitionParameters(windowsweep, type = type, ...)
  
  x.standardized <- (x - pp$mu.hat)/pp$s.hat
  if(plotme)
  {
    par(mfrow=c(1,3))
    qqnorm(x.standardized)
    qqline(x.standardized)
    hist(x.standardized, col="grey", breaks=50, freq=FALSE)
    curve(dnorm(x), col=2, lwd=2, add=TRUE)
    acf(x.standardized, na.action=na.pass)
    layout(1)
  }
  if(values) return(x.standardized)
}
