#' @title Visualize Distribution of Missing Values
#' 
#' @description Visualize the distribution of missing values within a time series. This is done
#' as a barplot, which is especially useful if the time series would otherwise be too large to be plotted.
#' 
#' @param x Numeric Vector (\code{\link{vector}}) or Time Series (\code{\link{ts}}) object containing NAs
#' 
#' @param colPoints Color of the points for each observation
#' @param colBackgroundMV Color for the background for the NA sequences
#' @param main Main label for the plot
#' @param xlab Label for x axis of the plot
#' @param ylab Label for y axis of plot
#' @param pch Plotting 'character', i.e., symbol to use.
#' @param cexPoints character (or symbol) expansion: a numerical vector.
#' @param col Color for the lines.
#' @param ... Additional graphical parameters that can be passed through to plot 
#' 
#' @details This function visualizes the distribution of missing values within a time series.
#' In comparison to the \code{\link[imputeTS]{plotNA.distribution}} function this is not done by plotting
#' each observation of the time series seperatly. Instead observations for time intervals are represented as bars.
#' For these intervals information about the amount of missing values are shown. This has the advantage, that also
#' for large time series a plot which is easy to overview can be created.
#' 
#' 
#' @author Steffen Moritz
#' 
#' @seealso \code{\link[imputeTS]{plotNA.distributionBar}},
#'  \code{\link[imputeTS]{plotNA.gapsize}}, \code{\link[imputeTS]{plotNA.imputations}}
#' 
#' @examples
#' #Prerequisite: Load a time series with missing values
#' x <- tsAirgap
#' 
#' #Example 1: Visualize the missing values in this time series
#' plotNA.distribution(x)
#' 
#' @importFrom graphics par plot points barplot
#' @export

plotNA.distribution <- function(x, colPoints = "steelblue", colBackgroundMV = "indianred2",
                                main="Distribution of NAs", xlab = "Time", ylab="Value", pch = 20, cexPoints = 0.8, col ="black", ... ) {
  
  data <- x
  
  #Check for wrong input 
  data <- precheck(data)

  id.na <- which(is.na(data))
  barplotData <- rep(NA,length(data))
  
  
  #make sure the end of the bar can not be seen
  barplotData[id.na] <- max(data, na.rm =T )*100
  
  ##Plot the red in background for unknown values
  barplot(barplotData, col = colBackgroundMV,xaxt = "n", yaxt = "n",   xlab = "", ylab = "", border = colBackgroundMV)
  
  ## Plot the line diagram of known values
  par(new=TRUE)
  plot(data, main =main, type = "l", xlab = xlab, ylab = ylab, col = col,... )
  points(data, pch= pch , cex = cexPoints, col = colPoints)
  
}




