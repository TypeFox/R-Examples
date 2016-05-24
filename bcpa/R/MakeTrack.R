#' Make Track
#'
#' Simple convenience function for creating a \code{track} class object from X-Y-Time movement data. A \code{track} class object can be conveniently plotted and analyzed within \code{bcpa}.
#' 
#' @param X vector of X locations
#' @param Y vector of Y locations
#' @param Time vector of time (can be POSIX)
#' @return a \code{track} class data frame, with three columns: \code{X}, \code{Y} and \code{Time}.
#' @seealso plot.track
#' @examples 
#' X <- cumsum(arima.sim(n=100, model=list(ar=0.8)))
#' Y <- cumsum(arima.sim(n=100, model=list(ar=0.8)))
#' Time <- 1:100
#' mytrack <- MakeTrack(X,Y,Time)
#' plot(mytrack)

MakeTrack <- function(X, Y, Time)
{
  MyTrack <- data.frame(Time = Time, X = X, Y = Y)
  class(MyTrack) <- c("track", "data.frame")
  return(MyTrack)
}