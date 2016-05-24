#' @export
#' @title Internal structure used by rstacks, rdeques, and rpqueues
#' @description
#' For use by rstacks and rdeques. Simply an environment with no parent,
#' reference for the data and the next node.
#' @param data data to reference with this node.
#' @return an environment.
#' 
rstacknode <- function(data) {
  newnode <- new.env(parent = emptyenv())
  newnode$data <- data
  newnode$nextnode <- NULL
  class(newnode) <- "rstacknode"
  return(newnode)
}