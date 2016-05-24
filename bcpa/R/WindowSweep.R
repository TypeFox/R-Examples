#' Perform window sweep for BCPA
#'
#' This is the workhorse function of the BCPA.  It performs a sweep of the time series, searching for most significant change points and identifying the parsimonious model according to an adjusted BIC.
#'
#' @param data the data to be analyzed.  Most typically, output of the \code{\link{GetVT}} function containing step lengths, absolute and turning angles, etc.
#' @param variable a character string representing the response to apply the BCPA to.  For example \code{"V*cos(theta)"}  for persistence velocity, or \code{"log(V)"} for log of velocity. 
#' @param windowsize integer size of the analysis window as a number of data points (not time units).  Should probably be no smaller than 20. 
#' @param windowstep integer step size of analysis. Values greater than 1 speed the analysis up. 
#' @param K senstivity parameter for the adjusted BIC.  Smaller values make for a less sensitive model selection, i.e. more likely that the null model of no significant changes will be selected.
#' @param tau a logical indicating whether the autocorrelation "rho" or the characteristic time "tau" should be estimated. 
#' @param range a number between 0 and 1 that determines the extent of each window that is scanned for changepoints.  I.e., if the window is 100 datapoints long, at the default \code{range=0.6}, changepoints will be scanned between 20 and 80. 
#' @param progress logical - whether or not to output a progress bar as the analysis is being performed. 
#' @param plotme logical - whether or not to plot the analysis as it is happening.  This slows the analysis considerably, especially in non-dynamic graphic environments like RStudio.   
#' @param ... additional parameters to be passed to the \code{\link{PartitionParameters}} function.
#' 
#' @return an object of class \code{windowsweep}, which is a list containing: 
#' \item{ws}{a data frame containing the change point, selected model, and parameter estimates} 
#' \item{x}{the response variable} 
#' \item{t}{the time variable - in the units specified in the data object} 
#' \item{t.POSIX}{the time variable as a POSIX objects (contianing Y-M-D H:S)} 
#' \item{windowsize}{the window size}
#' \item{windowstep}{the window step size}
#' 
#' @seealso  for internal functions: \code{\link{GetModels}}, \code{\link{GetBestBreak}}, \code{\link{GetDoubleL}}; for summarizing output: \code{\link{ChangePointSummary}}; for plotting output: \code{\link{plot.bcpa}}
#' @author Eliezer Gurarie
#' @examples
#' data(Simp)
#' plot(Simp)
#' Simp.VT <- GetVT(Simp)
#' Simp.ws <- WindowSweep(Simp.VT, "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
#' plot(Simp.ws, threshold=7)
#' plot(Simp.ws, type="flat", clusterwidth=3)
#' PathPlot(Simp, Simp.ws)
#' PathPlot(Simp, Simp.ws, type="flat")
#' DiagPlot(Simp.ws)

WindowSweep <- function (data, variable, windowsize = 50, windowstep = 1, K = 2, tau=TRUE, range=0.6, progress = TRUE, plotme=FALSE, ...) 
{
  x <- eval(parse(text = variable), data)
  t <- data$T.mid
  
  low <- seq(1, (length(t) - windowsize), windowstep)
  hi <- low + windowsize
  if(progress)
    pb <- txtProgressBar(min = 0, max = length(low), style = 3)
  
  for (i in 1:length(low)) {
    myx <- x[low[i]:hi[i]]
    myt <- t[low[i]:hi[i]]
    myestimate <- GetBestBreak(myx, myt, range, tau = tau)
    
    breakpoint <- myestimate[1]
    tbreak <- myestimate[2]
    allmodels <- GetModels(myx, myt, breakpoint, K, tau)
    
    # remember, column 3 is "bic"
    mymodel <- allmodels[allmodels[,3] == min(allmodels[,3]),]
    mymodel <- c(mymodel, Break = tbreak)
    
    if(i == 1)
      estimates <- mymodel
    else
      estimates <- rbind(estimates, mymodel)
    
    if(plotme)
    {
      plot.ts(t, x, type = "l", col = "grey")
      lines(t, x, type = "l")
      lines(myt, myx, col = "green")
      abline(v = tbreak)
    }
    
    # create progress bar
    if(progress & i %% 10 == 0) setTxtProgressBar(pb, i)
  }
  if(progress) close(pb)
  windowsweep <- list(ws = data.frame(estimates, row.names=1:nrow(estimates)), 
                      x=x, t=t, t.POSIX = data$T.POSIX, windowsize=windowsize, windowstep=windowstep)
  
  windowsweep$pp.smooth <- PartitionParameters(windowsweep, type="smooth", ...) 
  class(windowsweep) <- "bcpa"
  return(windowsweep)
}