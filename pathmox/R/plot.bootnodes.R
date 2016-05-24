#' @title Plot bootstrap results of terminal nodes
#' 
#' @description
#' Plot method for objects of class \code{"bootnodes"}
#' 
#' The function \code{plot.bootnodes} displays the value of the selected path
#' coefficient for the root node and the terminal nodes of PATHMOX and TECHMOX
#' trees. \cr In addition, the value of the mean bootstrap as well as the values
#' of the confidence interval are also shown.
#' 
#' @param x object of class \code{"bootnodes"}
#' @param pc integer indicating the number of path coefficient to be plotted
#' @param \dots Further arguments are ignored
#' @method plot bootnodes
#' @S3method plot bootnodes
plot.bootnodes <- function(x, pc = 1, ...)
{
  # ===================== Inputs =======================
  # PC = matrix of Path Coefficients
  # PMB = matrix of Mean values for Path Coeffs in Bootstrapping
  # PP05 = matrix of percentiles 0.05
  # PP95 = matrix of percentiles 0.95
  
  PC <- x$PC
  PMB <- x$PMB
  PP05 <- x$PP05
  PP95 <- x$PP95
  
  # Parameters
  min.y = min(PP05[pc,])# lower limit for axis 'y'
  max.y = max(PP95[pc,])# upper limit for axis 'y'
  
  # PLOT
  # dev.new()
  plot(1:ncol(PC), PC[pc,],  
       xlab = "Nodes", ylab=rownames(PC)[pc], 
       ylim = c(min.y - 0.1, max.y + 0.1), 
       main = c("Bootstrap Intervals for path coefficient", 
                rownames(PC)[pc], sep=" "),
       cex.main = 1, xaxt = "n", type = "n")
  abline(h=0, lty=2, col="grey")
  mtext(colnames(PC), side=1, at=1:ncol(PC), cex=.8)
  # plot of original path coefficients
  points(1:ncol(PC), PC[pc,], pch=20, col="blue")
  # plot of bootstrap path.coeffs
  points(1:ncol(PC), PMB[pc,], pch=21, cex=1.1, col="red")
  # plot of percentile 5%
  points(1:ncol(PC), PP05[pc,], pch="_", col="red")
  # plot of percentile 95%
  points(1:ncol(PC), PP95[pc,], pch="_", col="red")
  # plot of dotted line
  arrows(1:ncol(PC), PP05[pc,], 1:ncol(PC), PP95[pc,], 
         lty=3, col="red", code=0)
  legend(ncol(PC)-1, max.y + .1, c("Original", "Bootstrap"), 
         bty = "n", cex = 0.8, pch = c(19, 21), col = c("blue","red"))
}
