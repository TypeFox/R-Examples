
nodeNames  <- function(x) UseMethod("nodeNames")

nodeStates <- function(x, nodes=nodeNames(x)) UseMethod("nodeStates")

nodeNames.grain  <- function(x)
  x$universe$nodes

nodeStates.grain <- function(x, nodes=nodeNames(x)){
  x$universe$levels[nodes]
}
