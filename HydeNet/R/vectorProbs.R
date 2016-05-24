#' @name vectorProbs
#' @export vectorProbs
#' 
#' @title Convert a vector to JAGS Probabilities
#' @description Probability vectors can be passed manually to the model, but
#'   they must be formatted in code appropriate to JAGS.  \code{vectorProbs}
#'   will convert a vector of counts or weights to probabilities and format
#'   it into JAGS code.  
#'   
#' @param p a vector of counts, weights, or probabilities.
#' @param node the node for the parameters.  this is converted to a character
#'   string.  It is important that this be given accurately or it will not 
#'   match with the code written by \code{writeNetworkModel}.
#' @param normalize A logical indicating if the weights in \code{p} should
#'   be normalized (each value is taken as a proportion of the sum).
#'   
#' @author Jarrod Dalton and Benjamin Nutter
#' 
#' @examples
#' vectorProbs(c(1, 2, 3), "wells")
#'

vectorProbs <- function(p, node, normalize=TRUE){
  if (missing(node)) stop("Please provide the node name (note: spelling is important)")
  node <- as.character(substitute(node))
  if (normalize) p <- p/sum(p)
  paste0("pi.", node, "[", 1:length(p), "] <- ", p, collapse="; ")
}
