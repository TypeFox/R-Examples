#' Sankey / ribbon diagrams
#' 
#' Sankey diagrams are a type of flow diagrams, in which the width of the
#' arrows is proportional to the quantity they illustrate.
#' Riverplot allows the creation, in R, of a basic type of Sankey diagrams.
#'
#' First, you need to create a specific riverplot object that can be
#' directly plotted. (Use \code{\link{riverplot.example}} to generate an example object). 
#'
#' The simplest way is to create a graph-like representation of you diagram as
#' a list of nodes; each item in the list is a list of partner nodes.
#' Furthermore, you need to know at which position (from left to right) each
#' node resides.  Please take a look at the example section in the \code{\link{makeRiver}}
#' function.
#'
#' Once you have created a riverplot object with one of the above methods
#' (or manually), you can plot it either with \code{plot(x)} or
#' \code{riverplot(x)} (see \code{\link{riverplot}} for details).
#' 
#' @section Mini-gallery:
#' 
#' Simple example from \code{\link{riverplot.example}} function:
#' \code{plot( riverplot.example() )}.
#' 
#' \figure{example.jpg}{Simple example}
#'
#' Recreation of the famous figure by Charles Minard (see
#' \code{\link{minard}} for details).
#'
#' \figure{minard.jpg}{Minard}
#'
#'
#' @title Sankey / ribbon diagrams
#' @name riverplot-package
#' @docType package
#' @author January Weiner <january.weiner@@gmail.com>
NULL
