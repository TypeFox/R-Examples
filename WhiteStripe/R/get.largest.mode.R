#' @title Grab largest peak
#'
#' @description This function grabs the largest peak of the histogram
#' @param x values of midpoints from \code{\link{hist}}
#' @param y values of counts from \code{\link{hist}}
#' @param verbose print diagnostic output
#' @param ... arguments to be passed to \code{\link{smooth_hist}}
#' @export
#' @return Value of x that is the largest peak
#' @examples 
#' data(t2.voi.hist)
#' system.time({
#' y = t2.voi.hist$counts
#' x = t2.voi.hist$mids
#' x = x[!is.na(y)];
#' y = y[!is.na(y)]
#' # 30 used for speed of example
#' nawm_peak = get.largest.mode(x, y, k=30)
#' plot(t2.voi.hist, border="red")
#' abline(v=nawm_peak)
#' })
get.largest.mode <- function(x, y, verbose = TRUE,
  ...) {
  
  #estimate derivative
  if (verbose){
    cat("Smoothing Histogram\n")
  }
  system.time({
    smooth1 = smooth_hist(x, y, ...)
  })
  #estimate derivative
  if (verbose){
    cat("Smoothing Derivative\n")
  }  
  dy = get.deriv.smooth.hist(
    x,
    coefs=smooth1$coefs,
    knots=smooth1$knots,
    deg=smooth1$deg,
    deriv.deg=1)
  which.zero.crossing = which(
    (dy[1:(length(x)-1)]>0) > (dy[2:(length(x))]>0)
    )
  largest.peak = x[which.zero.crossing[which.max(y[which.zero.crossing])]]
  return(largest.peak)	
}

#' @title Get Last Peak
#'
#' @description This function grabs the last peak or shoulder.
#' @param x values of midpoints from \code{\link{hist}}
#' @param y values of counts from \code{\link{hist}}
#' @param rare.prop Proportion used to remove rare intensity tail
#' @param verbose print diagnostic output
#' @param remove.tail Remove rare intensity tail
#' @param ... arguments to be passed to \code{\link{smooth_hist}}
#' @export
#' @return Value of x that is the last peak
#' @examples
#' data(t1.voi.hist)
#' system.time({
#' y = t1.voi.hist$counts
#' x = t1.voi.hist$mids
#' x = x[!is.na(y)];
#' y = y[!is.na(y)]
#' # 20 used for speed of example
#' nawm_peak = get.last.mode(x, y, k=20)
#' plot(t1.voi.hist, border="red")
#' abline(v=nawm_peak)
#' })
#'  
get.last.mode = function(x,y, 
  rare.prop=1/5, verbose=TRUE, remove.tail = TRUE, ...) {
  
  
  #Remove rare intensity tail
  if (remove.tail){
    which.rare <- y < (rare.prop*max(y))
    y = y[!which.rare]
    x = x[!which.rare]
  }
  
  if (verbose){
    cat("Smoothing Histogram\n")
  }  
  #estimate derivative
  system.time({
    smooth1 = smooth_hist(x, y, ...)
  })
  #estimate derivative
  if (verbose){
    cat("Smoothing Derivative\n")
  }    
  dy<-get.deriv.smooth.hist(x, 
    coefs=smooth1$coefs,
    knots=smooth1$knots,
    deg=smooth1$deg,
    deriv.deg=1)
  which.zero.crossing<-which(
    (dy[1:(length(x)-1)]>0) > (dy[2:(length(x))]>0)
  )
  last.peak<-max(x[which.zero.crossing])
  #if (last.peak<median(t1.voi[t1.voi>mean(t1.voi)])) {
  # biggest.shoulder<-(x[x>median(t1.voi[t1.voi>mean(t1.voi)])])[which.min(abs(dy[x>median(t1.voi[t1.voi>mean(t1.voi)])]))]
  # return(biggest.shoulder)
  #} else {
  # return(last.peak)
  #}
  return(last.peak) 
}


#' @title Get First Peak
#'
#' @description This function grabs the first peak or shoulder.
#' @param x values of midpoints from \code{\link{hist}}
#' @param y values of counts from \code{\link{hist}}
#' @param rare.prop Proportion used to remove rare intensity tail
#' @param verbose print diagnostic output
#' @param remove.tail Remove rare intensity tail
#' @param ... arguments to be passed to \code{\link{smooth_hist}}
#' @export
#' @return Value of x that is the first peak
#' @examples
#' data(t1.voi.hist)
#' system.time({
#' y = t1.voi.hist$counts
#' x = t1.voi.hist$mids
#' x = x[!is.na(y)];
#' y = y[!is.na(y)]
#' # 20 used for speed of example
#' nawm_peak = get.first.mode(x, y, k=20)
#' plot(t1.voi.hist, border="red")
#' abline(v=nawm_peak)
#' })
#'  
get.first.mode = function(x,y, 
                         rare.prop=1/5, verbose=TRUE, remove.tail = TRUE, ...) {
  
  
  #Remove rare intensity tail
  if (remove.tail){
    which.rare <- y < (rare.prop*max(y))
    y = y[!which.rare]
    x = x[!which.rare]
  }
  
  if (verbose){
    cat("Smoothing Histogram\n")
  }  
  #estimate derivative
  system.time({
    smooth1 = smooth_hist(x, y, ...)
  })
  #estimate derivative
  if (verbose){
    cat("Smoothing Derivative\n")
  }    
  dy<-get.deriv.smooth.hist(x, 
                            coefs=smooth1$coefs,
                            knots=smooth1$knots,
                            deg=smooth1$deg,
                            deriv.deg=1)
  which.zero.crossing<-which(
    (dy[1:(length(x)-1)]>0) > (dy[2:(length(x))]>0)
  )
  last.peak<-min(x[which.zero.crossing])
  #if (last.peak<median(t1.voi[t1.voi>mean(t1.voi)])) {
  # biggest.shoulder<-(x[x>median(t1.voi[t1.voi>mean(t1.voi)])])[which.min(abs(dy[x>median(t1.voi[t1.voi>mean(t1.voi)])]))]
  # return(biggest.shoulder)
  #} else {
  # return(last.peak)
  #}
  return(last.peak) 
}
