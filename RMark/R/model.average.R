#' Compute model averaged estimates
#' 
#' A generic function to compute model averaged estimates and their standard
#' errors or variance-covariance matrix.
#' 
#' 
#' @aliases model.average model.average.default
#' @param x is either a list with a prescribed structure as defined in
#' \code{\link{model.average.list}} or a \code{marklist} as described in
#' \code{\link{model.average.marklist}}
#' @param ... additional arguments passed to specific functions
#' @return The structure of the returned value depends on which function is
#' called.
#' @author Jeff Laake
#' @export
#' @seealso
#' \code{\link{model.average.marklist}},\code{\link{model.average.list}}
#' @keywords utility
model.average<- function(x,...)
{
# A generic function that currently has definitions for default (model.average.default) and
# class marklist (model.average.marklist) 
UseMethod("model.average")
}
