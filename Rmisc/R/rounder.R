#' Round to Increment
#' 
#' Rounds a value to nearest increment
#' 
#' @param x The value to be rounded
#' @param inc The increment to round to
#' @param fun The rounding function. Valid options are 'floor', 'round' and 'ceiling'.
#' 
#' @return an object of class \code{numeric}
#' 
#' @export
#' 
#' @examples
#' rounder(.92, .05)
#' rounder(.93, .05)
#' rounder(.93, .05, "floor")
#' rounder(.93, .05, "ceiling")
#' 
rounder <- function(x, inc, fun="round")
{
  allowed <- c("floor","round","ceiling")
  if (fun %in% allowed)
    do.call(fun, list(x/inc)) * inc
  else
    stop("Invalid fun specified")
}