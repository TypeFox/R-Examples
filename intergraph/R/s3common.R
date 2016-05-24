
#============================================================================ 
# common functions

igVcount <- function(x, ...) UseMethod("igVcount")

igVcount.igraph <- function(x, ...)
  igraph::vcount(x)

igVcount.network <- function(x, ...)
  network::network.size(x)

igEcount <- function(x, ...) UseMethod("igEcount")

igEcount.igraph <- function(x, ...)
  igraph::ecount(x)

igEcount.network <- function(x, ...)
  network::network.edgecount(x)

igDirected <- function(x) UseMethod("igDirected")

igDirected.igraph <- function(x)
{
  igraph::is.directed(x)
}

igDirected.network <- function(x)
{
  network::is.directed(x)
}


