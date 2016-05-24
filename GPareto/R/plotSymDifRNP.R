##' Plot the symmetrical difference between two 
##' Random Non-Dominated Point (RNP) sets.
##' 
##' @title Symmetrical difference of RNP sets
##' @param set1,set2 RNP sets considered,
##' @param xlim,ylim numeric vectors of length 2, giving the \code{x} and \code{y} coordinates ranges for plotting,
##' @param fill optional color of the symmetric difference area,
##' @param add logical; if \code{TRUE} add to an already existing plot; if \code{FALSE} (default) start a new plot taking \code{xlim, ylim} as limits.
##' @param ... additional parameters for the \code{\link[graphics]{plot}} and \code{\link[graphics]{polygon}} graphic functions
##' @export
##' @examples
##' #------------------------------------------------------------
##' # Simple example
##' #------------------------------------------------------------
##' set1 <- rbind(c(0.2, 0.35, 0.5, 0.8),
##'               c(0.8, 0.6, 0.55, 0.3))
##' 
##' set2 <- rbind(c(0.3, 0.4),
##'               c(0.7, 0.4))
##' 
##' plotSymDifRNP(set1, set2, xlim = c(0, 1), ylim = c(0, 1), fill = "grey")
##' points(t(set1), col = "red", pch = 20)
##' points(t(set2), col = "blue", pch = 20)


plotSymDifRNP <- function(set1, set2, xlim, ylim, fill = "black", add = "FALSE", ...){
  
  #plot(t(cbind(set1,set2)),pch='.')
  lo <- nondominated_points(cbind(set1,set2))
  #points(t(lo),col="red")
  
  i1 <- length(set1[1,])
  i2 <- length(set2[1,])
  
  set1bis <- set1[,order(set1[1,])]
  set2bis <- set2[,order(set2[1,])]
  
  up <- -nondominated_points(-cbind(set1,set2,rbind(set1bis[1,2:i1],set1bis[2,1:(i1-1)]),
                                    rbind(set2bis[1,2:i2],set2bis[2,1:(i2-1)]),
                                    c(max(set1bis[1,1],set2bis[1,1]),ylim[2]),
                                    c(xlim[2],max(set1bis[2,i1],set2bis[2,i2]))
  )
  )
  #points(t(up),col="blue")
 
  #   plot(t(up))
  #   points(-t(lo),col="blue",pch=20)

  lo = lo[,order(lo[1,])]
  up = up[,order(up[1,])]
  
  if(add==FALSE){
    plot(NA,xlim=xlim,ylim=ylim, xlab = "", ylab = "", ...)
  }
  
  polygon(c(xlim[1],xlim[2],xlim[2],xlim[1]),c(ylim[1],ylim[1],ylim[2],ylim[2]),
          col = fill, border = NA, ...) 
  
  i1 <- length(lo[1,])
  i2 <- length(up[1,])
  
  # Add points to respect the dominance relationship
  lo1 <- cbind(lo,lo[,2:i1])
  
  lo1[,which(1:(2*i1)%%2 == 1)] = lo
  lo1[,which(1:(2*i1-1)%%2 == 0)] = rbind(lo[1,2:i1],lo[2,1:(i1-1)])
  
  
  
  polygon(t(cbind(lo1,c(xlim[2],lo[2,i1]),c(xlim[2],ylim[1]),c(xlim[1],ylim[1]),
                  c(xlim[1],ylim[2]),c(lo[1,1],ylim[2]))),col = "white", border = NA, ...)
  
  up1 <- cbind(up,up[,2:i2])
  
  up1[,which(1:(2*i2)%%2 == 1)] = up
  up1[,which(1:(2*i2-1)%%2 == 0)] = rbind(up[1,1:(i2-1)],up[2,2:i2])
  
  
  polygon(t(cbind(up1, c(xlim[2],up[2,i2]), c(xlim[2],ylim[2]),
                  c(up[1,1],ylim[2]))), col = "white", border = NA, ...)
  
  #return(list(lo=lo,up=up))
}