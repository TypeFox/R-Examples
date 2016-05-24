#' @export plotThin
#' @title Plot diagnosis for results of thin function 
#' 
#' @description
#' Three plots (selected by \code{which}) are currently available: 
#' a plot of the number of repetitions versus the number of maximum records retained
#' at each repetition ([1] observed values; [2] log transformed) and 
#' a histogram of the maximun records retained [3].
#' 
#' @param thinned A list of data.frames returned by \code{\link{thin}} function.
#' @param which if a subset of the plots is required, specify a subset of the numbers 1:3.
#' @param ask logical; if \code{TRUE}, the user is asked before each plot, see par(ask=.).
#' @param ... other parameters to be passed through to plotting functions.
#'
#' @seealso \code{\link{thin.algorithm}}
#' @seealso \code{\link{thin}}



plotThin <- function(thinned, which=c(1:3), 
                      ask=prod(par("mfcol")) < length(which)
                      && dev.interactive(), ...){
  
  par(ask=ask)
  
  ## Repetition number
  reps <- length(thinned)
  
  ## Look at the number of locs kept in each thinned dataset
  ## by determining the number of rows in each returned data.frame
  lat.long.thin.count <- unlist(lapply(thinned, nrow ))
  
  ## Create a vector of cummulative maximum records at each 
  ## repetition number
  cummax.lat.long.thin.count <- cummax(lat.long.thin.count)
  
  ## Plot the number of repetitions versus the number 
  ## of maximum records retained at each repetition
 if(any(1==which)){
  plot( (1:reps), cummax.lat.long.thin.count,
        xlab='Number Repetitions',
        ylab='Cummulative Maximum Records Retained',  
        xlim=c(0,reps), ...)
 }
  ## Make a log-log plot of the number of repetitions versus
  ## the number of maximum records retained
 if(any(2==which)){
 plot( log(1:reps), log(cummax.lat.long.thin.count),
        xlab='Log Number Repetitions',
        ylab='Log Cummulative Maximum Records Retained',
        #ylim=c(0,log(max(cummax.lat.long.thin.count))), 
        xlim=c(0,log(reps)), ...)
 }
  ## Plot a histogram of lat.long.thin.count
 if(any(3==which)){   
 hist(lat.long.thin.count,
      xlab='Maximum Records Retained',
      main="",
      ...)
 }
 
 par(ask=FALSE)
}