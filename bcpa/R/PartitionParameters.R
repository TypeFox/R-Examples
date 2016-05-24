#' Partition parameters 
#' 
#' Partitions - and, ultimately, estimates - all parameters of a BCPA, either as a rolling average (smooth BCPA), or as constant values within fixed change points (flat BCPA).
#' 
#' @param windowsweep a \code{windowsweep} object, i.e. the output of the \code{\link{WindowSweep}} function.
#' @param type type of estimate to output, whether "smooth" or "flat".
#' @param clusterwidth the temporal range within which changepoints are considered to be within the same cluster (for a "smooth" BCPA).
#' @return a data frame containing the three estimates: \code{mu.hat}, \code{s.hat}, \code{rho.hat}.
#' 
#' @author Eliezer Gurarie
#' @seealso used in \code{\link{ChangePointSummary}} and \code{\link{PhasePlot}}


PartitionParameters <-
function(windowsweep, type = c("smooth","flat")[1], clusterwidth = 1)
{
  if(type == "smooth" & !("pp.smooth" %in% names(windowsweep)))
  {
    x <- windowsweep$x
    t <- windowsweep$t 
    windowsize <- windowsweep$windowsize
    windowstep <- windowsweep$windowstep
    ws <- data.frame(windowsweep$ws)
    
    low <- seq(1, (length(t) - windowsize), windowstep)
    high <- low + windowsize
    
    n.col<-length(t)
    n.row<-dim(ws)[1]
    
    mu.M <- matrix(NA,n.row,n.col)
    s.M <- matrix(NA,n.row,n.col)
    rho.M <- matrix(NA,n.row,n.col)
    
    for(i in 1:n.row)
    {
      myws<-ws[i,]
      dts <- abs(t-myws$Break)
      tbreak <- match(min(dts),dts)
      
      max <- min(n.col,high[i])
      
      mu.M[i,low[i]:tbreak] <- myws$mu1
      mu.M[i,(tbreak+1):max] <- myws$mu2
      s.M[i,low[i]:tbreak] <- myws$s1
      s.M[i,(tbreak+1):max] <- myws$s2
      rho.M[i,low[i]:tbreak] <- myws$rho1
      rho.M[i,(tbreak+1):max] <- myws$rho2
    }
    
    adjust <- colSums(!is.na(mu.M))
    
    mu.hat<-colSums(mu.M,na.rm=1)/adjust
    s.hat<-colSums(s.M,na.rm=1)/adjust
    rho.hat<-colSums(rho.M,na.rm=1)/adjust
  } else
  {
    mu.hat <- windowsweep$pp$mu.hat
    s.hat <- windowsweep$pp$s.hat
    rho.hat <- windowsweep$pp$rho.hat
  }

  if(type == "flat")
  {
    x <- windowsweep$x 
    t <- windowsweep$t
    phases <- ChangePointSummary(windowsweep, clusterwidth)$phases             
 
    whichphase <- findInterval(t, phases$t0)
    
    mu.hat <- phases$mu.hat[whichphase]
    rho.hat <-  phases$rho.hat[whichphase]
    s.hat <-  phases$s.hat[whichphase]
  }

  return(data.frame(mu.hat,s.hat,rho.hat))
}
