#' Obtain summary of BCPA analysis
#'
#' Produces a summary of change points for a "flat" analysis, identifying phases (periods between change points) with estimated parameters, clustering neighboring ones according to a kernel density of the windowsweep breaks. 
#'
#' @param windowsweep a \code{windowsweep} object, i.e. the output of the \code{\link{WindowSweep}} function.
#' @param clusterwidth the temporal range within which change points are considered to be within the same cluster.  Corresponds to the bandwidth of the density of the break distribution.   
#' @param tau logical, whether to estimate the characteristic time tau (preferred) or not.  If FALSE, the autocorrelation parameter rho is calculated.  

#' @return a list containing two elements: 
#' \item{breaks}{a data frame containing columns:  \code{middle} - each change point, \code{size} - the number of windows that selected the change point, \code{modelmode} - the most frequently selected of the seven possible models (M0 - the null model - is excluded), and \code{middle.POSIX} - the mid-point as a POSIX time object.} 
#' \item{phases}{a data frame containing columns \code{mu.hat, s.hat, rho.hat} - the estimated mean, standard deviation, and time scale (or auto-correlation) within each phase (i.e. period between change points), \code{t0} - the beginning of the phase, \code{t1} - the end of the phase, and \code{interval} - the total duration of the phase.}
#' 
#' @author Eliezer Gurarie
#' @seealso \code{\link{WindowSweep}}
#' @examples
#' if(!exists("Simp.VT")){
#'  data(Simp)
#'  Simp.VT <- GetVT(Simp)}
#' if(!exists("Simp.ws"))
#'  Simp.ws <- WindowSweep(Simp.VT, "V*cos(Theta)", windowsize = 50, windowstep = 1, progress=TRUE)
#' # too many change points:
#' ChangePointSummary(Simp.ws)
#' # about the right number of change points:
#' ChangePointSummary(Simp.ws, clusterwidth=3)

ChangePointSummary <- function(windowsweep, clusterwidth=1, tau=TRUE)
{  
  ws <- subset(windowsweep$ws, windowsweep$ws$Model > 0)
  
  x <- windowsweep$x
  t <- windowsweep$t
  
  # Use density to obtain 
  breaks <- ws$Break
  break.density <- density(breaks, bw=clusterwidth)
  cps <- break.density$x[which(diff(diff(break.density$y) > 0) == -1) + 1]
  BreakCluster <- alply(t(cps), 2, function(cp) which(breaks < cp+clusterwidth & breaks > cp-clusterwidth))
   
  # THIS IS A LITTLE NUTTY!  THERE *MUST* BE A BETTER WAY TO GET THE MOST FREQUENT INCIDENCE!
  getMode <- function(x) as.numeric(names(sort(table(x), decreasing=TRUE)[1]))
  BreakTable <- ldply(BreakCluster, function(x) if(length(x) > 0)	
    data.frame(middle = mean(ws$Break[x]), size = length(ws$Break[x]), modelmode = 	getMode(ws$Model[x])))
  
  # obtain POSIX time of change points
  t.POSIX <- windowsweep$t.POSIX
  t.breaks <- c(min(t)-1, BreakTable$middle, max(t))  
  t.cut <- cut(t, t.breaks)
  
  index.breaks <- which(diff(as.numeric(t.cut)) == 1)
  break.POSIX1 <- t.POSIX[index.breaks]
  break.POSIX2 <- t.POSIX[index.breaks+1]
  if("POSIXt" %in% class(break.POSIX2))
    BreakTable$middle.POSIX <- break.POSIX1 + difftime(break.POSIX2 , break.POSIX1)/2  else BreakTable$middle.POSIX <- (break.POSIX1 + break.POSIX2)/2
  
  # Partitioning of periods
  phases <- ddply(data.frame(x,t), .(t.cut), 
                  function(xt) c(mu.hat = mean(xt$x), s.hat = sd(xt$x), rho.hat = GetRho2(xt$x, xt$t, tau=tau)[1]))
  phases <- data.frame(phases, t0 = t.breaks[-length(t.breaks)], 
                       t1 = t.breaks[-1],
                       interval = diff(t.breaks))
  
  return(list(breaks = BreakTable, phases = phases))
}


##########
# DEPRECATE BECAUSE DEPENDS ON INTERVALS
##########

# ChangePointSummary2 <- function(windowsweep, clusterwidth=1, tau=TRUE)
# {  
#   ws <- subset(windowsweep$ws, windowsweep$ws$Model > 0)
#   
#   x <- windowsweep$x
#   t <- windowsweep$t
#   
#   # Summary of breaks
#   
#   BreakCluster <- clusters(ws$Break, w=clusterwidth, which=TRUE)  
#   
#   # THIS IS A LITTLE NUTTY!  THERE *MUST* BE A BETTER WAY TO GET THE MOST FREQUENT INCIDENCE!
#   getMode <- function(x) as.numeric(names(sort(table(x), decreasing=TRUE)[1]))
#   
#   BreakTable <- ldply(BreakCluster, function(x) 	
#     data.frame(middle = mean(ws$Break[x]), size = length(ws$Break[x]), modelmode = 	getMode(ws$Model[x])))
#   
#   # obtain POSIX time of change points
#   t.POSIX <- windowsweep$t.POSIX
#   t.breaks <- c(min(t)-1, BreakTable$middle, max(t))  
#   t.cut <- cut(t, t.breaks)
#   
#   index.breaks <- which(diff(as.numeric(t.cut)) == 1)
#   break.POSIX1 <- t.POSIX[index.breaks]
#   break.POSIX2 <- t.POSIX[index.breaks+1]
#   if("POSIXt" %in% class(break.POSIX2))
#     BreakTable$middle.POSIX <- break.POSIX1 + difftime(break.POSIX2 , break.POSIX1)/2  else BreakTable$middle.POSIX <- (break.POSIX1 + break.POSIX2)/2
#   
#   # Partitioning of periods
#   phases <- ddply(data.frame(x,t), .(t.cut), 
#                   function(xt) c(mu.hat = mean(xt$x), s.hat = sd(xt$x), rho.hat = GetRho2(xt$x, xt$t, tau=tau)[1]))
#   phases <- data.frame(phases, t0 = t.breaks[-length(t.breaks)], 
#                        t1 = t.breaks[-1],
#                        interval = diff(t.breaks))
#   
#   return(list(breaks = BreakTable, phases = phases))
# }
