#' 'Flattening' a list-form crwPredict object into a data.frame
#' 
#' \dQuote{Flattens} a list form \code{\link{crwPredict}} object into a flat
#' data.frame.
#' 
#' 
#' @param predObj A crwPredict object
#' @return a \code{\link{data.frame}} version of a crwPredict list with columns
#' for the state standard errors
#' @author Devin S. Johnson
#' @seealso \code{\link{northernFurSeal}} for use example
#' @export
"as.flat" <- function(predObj)
{
    se.y <- sqrt(t(apply(predObj$V.hat.y, 3, diag)))
    se.x <- sqrt(t(apply(predObj$V.hat.x, 3, diag)))
    colnames(se.y) <- paste("se", names(predObj$alpha.hat.y), sep=".")
    colnames(se.x) <- paste("se", names(predObj$alpha.hat.x), sep=".")
    flat <- cbind(predObj$originalData, predObj$alpha.hat.y, se.y,
                  predObj$alpha.hat.x, se.x)
    if (!is.null(predObj$speed)) flat <- cbind(flat, predObj$speed)
    class(flat) <- c("crwPredict", "data.frame")
    attr(flat, "coord") <- attr(predObj, "coord")
    attr(flat, "random.drift") <- attr(predObj, "random.drift")
    attr(flat, "stop.model") <- attr(predObj, "stop.model")
    attr(flat, "polar.coord") <- attr(predObj, "polar.coord")
    attr(flat, "Time.name") <- attr(predObj, "Time.name")
    attr(flat, "flat") <- TRUE
    return(flat)
}

#' 'Flattening' a list-form crwPredict object into a data.frame
#' 
#' \dQuote{Flattens} a list form \code{\link{crwPredict}} object into a flat
#' data.frame.
#' 
#' 
#' @param predObj A crwPredict object
#' @return a \code{\link{data.frame}} version of a crwPredict list with columns
#' for the state standard errors
#' @author Devin S. Johnson
#' @seealso \code{\link{northernFurSeal}} for use example
#' @export
"flatten" <- function(predObj)
{
  se <- sqrt(t(apply(predObj$V.hat, 3, diag)))
  colnames(se) <- paste("se", names(predObj$alpha.hat), sep=".")
  flat <- cbind(predObj$originalData, predObj$alpha.hat, se)
  if (!is.null(predObj$speed)) flat <- cbind(flat, speed=predObj$speed)
  class(flat) <- c("crwPredict", "data.frame")
  attr(flat, "coord") <- attr(predObj, "coord")
  attr(flat, "random.drift") <- attr(predObj, "random.drift")
  attr(flat, "Time.name") <- attr(predObj, "Time.name")
  attr(flat, "flat") <- TRUE
  return(flat)
}
