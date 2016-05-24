#' @family Wrapper Generators
#' 
#' @title Create a boostr compatible wrapper for a reweighter.
#' @description Use provided metadata on a given reweighter to create a boostr
#' compatible wrapper.
#' 
#' @param reweighter a function which satisfies the abstract definition of a
#' reweighter (see description below).
#' @param reweighterInputPreds a string indicating the name of the argument 
#' \code{reweighter} uses to represent the input predictions.
#' @param reweighterInputResponse a string indicating the name of the argument 
#' \code{reweighter} uses to represent the true responses for the input predictions.
#' @param reweighterInputWts a string indicating the name of the argument 
#' \code{reweighter} uses to represent the input weights.
#' @param reweighterOutputWts a string indicating the name of the entry in
#' \code{reweighter}'s output that represents the output weights.
#' @param .verbose a boolean indicating if warnings should be displayed or not.
#' 
#' @template reweighters
#' 
#' @return A function (wrapper around \code{reweighter}) which is a '\code{reweighter}'
#' object. The wrapper's signature is \code{(prediction, response, weights, ...)}
#' and its output is a list that names the cell containing its weights
#' '\code{weight}'. Hence, the wrapper is a boostr compatible reweighter.
#'



wrapReweighter <- function(reweighter,
                        reweighterInputPreds="prediction",
                        reweighterInputResponse="response",
                        reweighterInputWts="weights",
                        reweighterOutputWts="weights",
                        .verbose=FALSE) {
  # extend signature of reweighter
  reweighter <- addDots(reweighter, .verbose)
  f <- function(prediction, response, weights, ...) {
    reweighterArgs <- list(prediction=prediction, response=response,
                        weights=weights)
    
    names(reweighterArgs)[1:3] <- c(reweighterInputPreds, reweighterInputResponse,
                                 reweighterInputWts)
    
    output <- do.call(reweighter, c(reweighterArgs, list(...)))
    
    if (!inherits(output, 'list') && .verbose) {
      warning("reweighter output is not a list -- this may cause problems")
    }
    
    names(output)[names(output) == reweighterOutputWts] <- "weights"
    
    output
  }
  class(f) <- c("reweighter", class(f))
  f
}

