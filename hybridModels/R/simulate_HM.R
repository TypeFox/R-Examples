#' It runs the chosen hybrid model.
#' 
#' @description  \code{simHM} is generic function that calls a method to run the
#'               simulation base on object's class
#' 
#' @param x of a specific class of model.
#' 
#' @inheritParams hybridModel
#' 
#' @return A \code{\link{data.frame}} with the number of individuals through
#'         time per node, per state and per simulation.
#'
#' @references  .
#' @seealso \link{GillespieSSA}.
#' @export
#' @import foreach
simHM <- function(x, network, sim.number, num.cores = 'max', fill.time) UseMethod("simHM")