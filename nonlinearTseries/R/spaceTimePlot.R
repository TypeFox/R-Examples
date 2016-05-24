################################################################################
#' Space Time plot
#' @description
#' The space time separation is a broadly-used method of detecting non-stationarity
#' and temporal correlations in the time series being analyzed. The space time
#' separation plot is also used to select a proper Theiler window by selecting 
#' a temporal separation enough to saturate the contour lines.
#' @details
#' Each contour line of the space time plot indicate the distance you have to go (y-axis) to find a given fraction of 
#' neighbour pairs, depending on their temporal separation (x-axis).
#' 
#' WARNING: The parameter \emph{number.time.steps} should be used with caution since this 
#' method performs heavy computations.
#' @param time.series The original time series being analyzed.
#' @param embedding.dim Integer denoting the dimension in which we shall embed the time series.
#' @param time.lag Integer denoting the number of time steps that will be use to construct the 
#' Takens' vectors.
#' @param takens Instead of specifying the \emph{time.series}, the \emph{embedding.dim} and the \emph{time.lag}, the user
#' may specify directly the Takens' vectors. 
#' @param max.radius Maximum neighbourhood radius in which the algorithm will look for finding neighbours. This
#' parameter may be used to avoid heavy computations. If the user does not specify a radius, the algorithm estimates it.
#' @param time.step Integer denoting the number of discrete steps between two calculations of the space time plot.
#' @param number.time.steps Integer denoting the number of temporal jumps in steps of \emph{time.step} in which we want to compute the space time plot.
#' @param numberPercentages Number of contour lines to be computed. Each contour line represent a concrete percentage of points (see Details).
#' @param do.plot Logical. If TRUE, the time space plot is shown.
#' @param ... Additional plotting parameters.
#' @return A \emph{timeSpacePlot} object that consist, essentially, of a matrix storing the values for each contour line.
#' Each row stores the value for a given percentage of points. Each column stores the value of the radius you have to go to
#'  find a given fraction of neighbour pairs (the rows), depending on their temporal separation (the colums). This matrix 
#'  can be accessed by using the \emph{contourlines} method.
#' @references  H. Kantz  and T. Schreiber: Nonlinear Time series Analysis (Cambridge university press)
#' @examples
#'  \dontrun{
#'  tak = buildTakens(sin(2*pi*0.005*(0:5000)),2,1)
#'  stp.test = spaceTimePlot(takens=tak,number.time.steps=400,do.plot=TRUE)
#'  }
#' @author Constantino A. Garcia
#' @rdname spaceTimePlot
#' @export spaceTimePlot
#' @exportClass spaceTimePlot
#' @useDynLib nonlinearTseries
spaceTimePlot=function(takens = NULL, 
                       time.series=NULL, embedding.dim=2, 
                       time.lag=1,max.radius=NULL, time.step=1,
                       number.time.steps = NULL, numberPercentages=10,
                       do.plot=TRUE,...){
  ### define
  kLengthRadiusVector = 1000
  kTimeStepsDefault = 500
  ###
  if(is.null(takens)){
    takens = buildTakens( time.series, embedding.dim = embedding.dim, time.lag = time.lag)  
  } 
  if(is.null(number.time.steps)){
    number.time.steps = min(length(takens)/time.step,kTimeStepsDefault)
  } 
    
  if (is.null(max.radius)){max.radius = (max(takens)-min(takens)+0.01)}
  # parameters for the C function
  numberTakens = nrow(takens)
  embedding.dim = ncol(takens)
  radius =  ( (1:kLengthRadiusVector)* max.radius )/kLengthRadiusVector
  stp = matrix(0,ncol=number.time.steps,nrow=numberPercentages)
  # call the C function
  sol = .C("spaceTimePlot",takens = as.double(takens), numberTakens = as.integer(numberTakens), 
                   embeddingD = as.integer(embedding.dim), eps = as.double(radius), leps = as.integer(kLengthRadiusVector),
                   numberTimeSteps =  as.integer(number.time.steps), timeStep = as.integer(time.step),
                   numberPercentages = as.integer( numberPercentages), spaceTimePlot = as.double(stp),
                   PACKAGE="nonlinearTseries")
  # format and return spaceTimePlot
  stp.matrix = matrix(sol$spaceTimePlot,ncol=number.time.steps,nrow=numberPercentages)
  # positions where the radius was not enough to compute the propper percentage 
  positions = which(stp.matrix==0,arr.ind=TRUE)
  if (length(positions) == 0){
    warning("The maximum radius was not enough to find enough neighbours for all the percentages\n")
    stp.matrix[positions]=NA  
  }
  
  dimnames(stp.matrix)  = list(percentagePoints=100*(1:numberPercentages)/numberPercentages,number.time.steps=1:number.time.steps)
  stp = list(stp.matrix=stp.matrix,time.step=time.step,time.axis = 1:number.time.steps)
  class(stp) = "spaceTimePlot"
  stp = propagateTakensAttr(stp,takens)
  
  #plot if necessary
  if (do.plot){
    plot(stp,...)
  }
  
  return (stp)
  
}


#' Obtain the contour lines of the space time plot.
#' @param x A \emph{spaceTimePlot} object.
#' @return Returns a matrix representing the contour lines of the 
#' space time plot.
#' @seealso \code{\link{spaceTimePlot}}
#' @export contourLines
contourLines = function(x){
  UseMethod("contourLines")
}


#' @return The \emph{contourLines} function returns the contour lines of the 
#' space time plot.
#' @param x A \emph{spaceTimePlot} object.
#' @rdname spaceTimePlot
#' @method contourLines spaceTimePlot
#' @export
contourLines.spaceTimePlot = function(x){
  return (x$stp.matrix)
}
  

#' @rdname spaceTimePlot
#' @method plot spaceTimePlot
#' @param main A title for the plot.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param type Type of plot (see \code{\link[graphics]{plot}}).
#' @param ylim Numeric vector of length 2, giving the y coordinates range.
#' @param col Vector of colors for each of the percentages of the plot.
#' @param pch Vector of symbols for each of the percentages of the plot.
#' @param add.legend add a legend to the plot?
#' @export
plot.spaceTimePlot = function(x, main = "Space time separation plot",xlab=NULL,
                              ylab=NULL,type = "l", ylim=NULL,
                              col = NULL, pch = NULL,
                              add.legend=TRUE,
                              ...){
  # auxiliar parameters for plotting
  timeSteps=x$time.axis
  numberPercentages = nrow(x$stp.matrix)
  if (is.null(xlab)) xlab = paste("Number of time steps, in steps of ",
                                  x$time.step)
  if (is.null(ylab)) ylab = expression("Distance to reference point ("*epsilon*")")
  if (is.null(ylim)) ylim = ylim=c(0,1.1*max(x$stp.matrix))
  # obtain vector of graphical parameters if not specified
  col = vectorizePar(col,numberPercentages)
  pch = vectorizePar(pch,numberPercentages)
  
  if (add.legend){
    current.par =  par(no.readonly = TRUE)
    on.exit(par(current.par))
    layout(rbind(1,2), heights=c(8,2))
  }
  # plot first using plot
  plot(timeSteps,x$stp.matrix[1,],type=type,col=col[[1]],pch=pch[[1]]
       ,xlab=xlab, ylab = ylab, main =main, ylim=ylim,...)
  legend.text = paste(100*1/numberPercentages,"%")
  # plot using lines
  if (numberPercentages>1){
    for (i in 2:numberPercentages){
      lines(timeSteps,x$stp.matrix[i,],type=type, col=col[[i]],pch=pch[[i]],...)
      legend.text = c(legend.text,paste(100*i/numberPercentages,"%"))
    }
  }
  if(add.legend){
    par(mar=c(0, 0, 0, 0))
    # c(bottom, left, top, right)
    plot.new()
    legend("center","groups",ncol=numberPercentages/2,
           col=col,lty=rep(1,numberPercentages),
           lwd=rep(2.5,numberPercentages),bty="n",
           legend=legend.text,title="Percentage of neighbour pairs")

  }
  
}
