#' Plot set of pixels on grid
#' 
#' \code{plot} plots a grid of pixels based on an object from \code{pgrid} or \code{confreg}.
#' 
#' If a vector of pixel indices is supplied to \code{set}, then those pixels will be colored \code{col} by this function and the \code{type} argument has no effect.  On the other hand, if the \code{set} argument is of class \code{confreg}, then the function digs in to display either the \code{confidence} or \code{complement} set in the \code{confreg} object.  In that case, \code{type} is used to decide which set to display.
#' 
#' @param x An object returned from the \code{pgrid} function.
#' @param set A vector which contains the indices of the pixels/cells that should be plotted.  OR a \code{confreg} object from the \code{confreg} function.  See Details.
#' @param col The color of the plotted pixels.
#' @param add A logical value indicating whether the pixels should be added to an existing plot (\code{add = TRUE}) or should the pixels be plotted on a new plot (\code{add = FALSE}).
#' @param type The type of set of plot if \code{set} of of class \code{confreg}.  Th default is \code{"confidence"}, while the other option is \code{complement}, based on the components of the \code{confreg} object.
#' @param ... Additional arguments that will be passed to the plot function (assuming \code{add} \code{=} \code{FALSE}.)
#' 
#' @return This function does not return anything; it only creates a new plot or modifies an existing plot.
#' @author Joshua French
#' @export
#' @examples
#' library(SpatialTools)
#' 
#' # Example for exceedance regions
#' 
#' set.seed(10)
#' # Load data
#' data(sdata)
#' # Create prediction grid
#' pgrid <- create.pgrid(0, 1, 0, 1, nx = 26, ny = 26)
#' pcoords <- pgrid$pgrid
#' # Create design matrices
#' coords = cbind(sdata$x1, sdata$x2)
#' X <- cbind(1, coords)
#' Xp <- cbind(1, pcoords)
#' 
#' # Generate covariance matrices V, Vp, Vop using appropriate parameters for 
#' # observed data and responses to be predicted
#' spcov <- cov.sp(coords = coords, sp.type = "exponential", 
#'  sp.par = c(1, 1.5), error.var = 1/3, finescale.var = 0, pcoords = pcoords)
#' 
#' # Predict responses at pgrid locations
#' krige.obj <- krige.uk(y = as.vector(sdata$y), V = spcov$V, Vp = spcov$Vp, 
#'  Vop = spcov$Vop, X = X, Xp = Xp, nsim = 100, 
#'  Ve.diag = rep(1/3, length(sdata$y)) , method = "chol")
#'                 
#' # Simulate distribution of test statistic for different alternatives
#' statistic.sim.obj.less <- statistic.sim(krige.obj = krige.obj, level = 5, 
#'  alternative = "less")
#' statistic.sim.obj.greater <- statistic.sim(krige.obj = krige.obj, 
#'  level = 5, alternative = "greater")
#' # Construct null and rejection sets for two scenarios
#' n90 <- exceedance.ci(statistic.sim.obj.less, conf.level = .90, 
#'  type = "null")
#' r90 <- exceedance.ci(statistic.sim.obj.greater,conf.level = .90, 
#'  type = "rejection")       
#' # Plot results
#' plot(pgrid, n90, col="blue", add = FALSE, xlab = "x", ylab = "y")
#' plot(pgrid, r90, col="orange", add = TRUE)
#' legend("bottomleft", 
#'  legend = c("contains true exceedance region with 90 percent confidence", 
#'  "is contained in true exceedance region with 90 percent confidence"),
#'  col = c("blue", "orange"), lwd = 10)  


plot.pgrid <- function(x, set, col = "gray", add = FALSE, type = "confidence", ...)
{
  nx <- length(x$upx)
  ny <- length(x$upy)
  if(class(set) == "confreg")
  {
    type = match.arg(type, c("confidence", "complement"))
    conf <- comp <- rep(NA, nx*ny)
    conf2 <- comp2 <- rep(NA, nrow(x$pgrid))
    conf2[set$conf] <- 1
    comp2[set$comp] <- 1
    conf[x$p.in.grid] <- conf2
    comp[x$p.in.grid] <- comp2
    
    if(type == "confidence"){ setvec = conf }
    else{ setvec = comp } 
    image(x$upx, x$upy, matrix(setvec, nx, ny), col = col, add = add, ...)
  }
  else
  {
    setvec <- rep(0, nrow(x$pgrid))
    setvec[set] <- 1
    image(x$upx, x$upy, matrix(setvec, nx, ny), col = c(0, col), add = add, ...)
  }
}



