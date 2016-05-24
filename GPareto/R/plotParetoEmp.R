##' Plot the Pareto front with step functions.
##' 
##' @title Pareto front visualization
##' @param nondominatedPoints points considered to plot the Pareto front with segments, matrix with one point per row,
##' @param add optional boolean indicating whether a new graphic should be drawn,
##' @param max optional boolean indicating whether to display a Pareto front in a maximization context,
##' @param ... additional values to be passed to the \code{\link[graphics]{lines}} function.
##' @export
##' @examples
##' #------------------------------------------------------------
##' # Simple example
##' #------------------------------------------------------------
##' 
##' x <- c(0.2, 0.4, 0.6, 0.8)
##' y <- c(0.8, 0.7, 0.5, 0.1)
##' 
##' plot(x, y, col = "green", pch = 20) 
##' 
##' plotParetoEmp(cbind(x, y), col = "green")
##' ## Alternative
##' plotParetoEmp(cbind(x, y), col = "red", add = FALSE)
##' 
##' ## With maximization
##' 
##' plotParetoEmp(cbind(x, y), col = "blue", max = TRUE)

plotParetoEmp <- function(nondominatedPoints, add = TRUE, max = FALSE,...){
  if(class(nondominatedPoints) != "matrix"){
    cat("The nondominatedPoints argument should be a matrix \n")
  }else{
    if(add == FALSE){
      plot(nondominatedPoints, ...)
    }
    
    temp <- nondominatedPoints[order(nondominatedPoints[,1]), , drop = FALSE]
    
    if(max){
      # Limiting points
      if(add){
        lim_left <- par('usr')[1]
        lim_bottom <- par('usr')[3]
      }else{
        lim_left <- -abs(20*temp[nrow(temp),1])
        lim_bottom <- -abs(20*temp[1,2])
      }
      
      lines(c(lim_left, temp[1, 1]), c(temp[1, 2], temp[1, 2]), ...)
      
      lines(c(temp[dim(temp)[1], 1], temp[dim(temp)[1], 1]), c(temp[dim(temp)[1], 2], lim_bottom), ...)
      
      # Segments in between
      if(nrow(temp) > 1){
        for (i in 1:(nrow(temp) - 1)) {
          lines(c(temp[i, 1], temp[i , 1]), c(temp[i, 2], temp[i+ 1, 2]), ...)
          lines(c(temp[i, 1], temp[i + 1, 1]), c(temp[i+1, 2],  temp[i + 1, 2]), ...)
        }
      }
      
    }else{
      # Limiting points
      if(add){
        lim_right <- par('usr')[2]
        lim_up <- par('usr')[4]
      }else{
        lim_right <- abs(20*temp[nrow(temp),1])
        lim_up <- abs(20*temp[1,2])
      }
      
      lines(c(temp[1,1], temp[1,1]),c(lim_up, temp[1,2]),...) #first segment
      
      lines(c(temp[dim(temp)[1], 1], lim_right),
            c(temp[dim(temp)[1], 2], temp[dim(temp)[1], 2]),...) #last segment
      
      # Segments in between
      if(nrow(temp) > 1){
        for(i in 1:(nrow(temp)-1)){
          lines(c(temp[i,1], temp[i+1,1]),c(temp[i,2], temp[i,2]),...) #horizontal part
          lines(c(temp[i+1,1], temp[i+1,1]),c(temp[i,2], temp[i+1,2]),...) #vertical part
        }
      }
    }
  }
}
