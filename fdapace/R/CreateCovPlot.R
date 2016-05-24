#' Create the covariance surface plot based on the results from FPCA() or FPCder().
#'
#' This function will open a new device if not instructed otherwise.
#'
#' @param fpcaObj returned object from FPCA().
#' @param covPlotType a string specifying the type of covariance surface to be plotted:
#'                     'Smoothed': plot the smoothed cov surface 
#'                     'Fitted': plot the fitted cov surface
#' @param isInteractive an option for interactive plot:
#'                      TRUE: interactive plot; FALSE: printable plot
#' @param ... other arguments passed into persp3d, persp3D, plot3d or points3D for plotting options
#'
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' CreateCovPlot(res)
#' @export

CreateCovPlot = function(fpcaObj, covPlotType = 'Fitted', isInteractive = FALSE, ...){
  
  ## Check if plotting covariance surface for fitted covariance surface is proper
  if(covPlotType == 'Fitted'){
    no_opt = length(fpcaObj$lambda)
    if(length(no_opt) == 0){
      warning('Warning: Input is not a valid FPCA or FPCder output.')
      return()
    } else if(no_opt == 1){
      warning('Warning: Covariance surface is not available when only one principal component is used.')
      return()
    }
  }
  
  ## Check if rgl is installed
  if(isInteractive == TRUE && is.element('rgl', installed.packages()[,1]) == FALSE){
    isInteractive = FALSE
    warning("Interactive plot requires package 'rgl', isInteractive set to be FALSE!")
  }
  
  ## Define the variables to plot
  workGrid = fpcaObj$workGrid
  if(covPlotType == 'Smoothed'){
    covSurf = fpcaObj$smoothedCov # smoothed covariance matrix
  } else if(covPlotType == 'Fitted'){
    covSurf = fpcaObj$fittedCov # fitted covariance matrix
  } else {
    warning("Covariance plot type no recognised; using default ('Fitted').")
    covSurf = fpcaObj$fittedCov 
    covPlotType = 'Fitted'
  }
  
  ## Define the plotting arguments
  args1 <- list( xlab='s', ylab='t', col =  topo.colors(24), zlab = 'C(s,t)', lighting=TRUE)
  if( covPlotType == 'Fitted' ){
    args1$main = 'Fitted covariance surface';
  } else {
    args1$main = 'Smoothed covariance surface';
  }
  if (isInteractive){ # If is interactive make it blue instead of topological
    args1$col = 'blue'
  }
  inargs <- list(...)
  args1[names(inargs)] <- inargs 
  args2 = list (x = workGrid, y = workGrid, z = covSurf)
  
  ## Plot the thing
  if(isInteractive){ # Interactive Plot 
    do.call(rgl::persp3d, c(args2, args1)) 
  } else { # Static plot 
    do.call(plot3D::persp3D, c(args2, args1))    
  }  
}
