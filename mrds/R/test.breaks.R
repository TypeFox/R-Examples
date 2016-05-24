#' Test validity for histogram breaks(cutpoints)
#'
#' Determines whether user specified breaks for histograms are properly ordered
#' and match the left and right truncation.
#'
#'
#' @param breaks vector of cutpoints (breaks) for distance histogram
#' @param left left truncation value
#' @param width right truncation value; either radius of point count or
#'   half-width of transect
#' @return vector of breaks modified to be valid if necessary
#' @author Jeff Laake
test.breaks <- function(breaks,left,width){
  ### Define function make.bins
  ##make.bins <- function(xmat,bins){
  ##  indices <- as.numeric(cut(xmat$distance,bins,include.lowest=TRUE))
  ##  xmat$distbegin <- bins[indices]
  ##  xmat$distend <- bins[indices+1]
  ##  return(xmat)
  ##}

  # Make sure break points are in order
  if(any(breaks!=sort(breaks))){
    stop("Break points are out of order.")
  }

  # if any endpoint > width, issue warning and reset endpoints
  if(any(breaks>1.000001*width)){
    message(paste("Specified endpoints > ",width,"; values reset."))
    breaks <- c(breaks[breaks<width],width)
  }

  # if last endpoint does not include width reset last endpoint
  if(breaks[length(breaks)]<0.999999*width){
    message(paste("Last interval endpoint did not include", width,
                  ". It was reset."))
    breaks <- c(breaks,width)
  }

  # if any endpoint < left, issue warning and reset endpoints
  if(any(breaks<0.99999*left)){
    message(paste("Specified endpoints < ",left,"; values reset."))
    breaks <- c(left,breaks[breaks>left])
  }

  # if first endpoint does not include left reset last first point
  if(breaks[1]>1.00001*left){
    message(paste("First interval endpoint did not include,",left,
                  ". It was reset"))
     breaks <- c(left,breaks)
  }
  return(breaks)
}
