#' @import igraph Matrix lattice gdata XML cluster clusterSim
NULL

#' Constructing a conceptmaps object
#' 
#' \code{conceptmaps} creates an object that encompasses a set of concept maps. For actual implementations, see
#' \code{\link{conceptmaps.default}}, \code{\link{conceptmaps.list}}, \code{\link{conceptmaps.matrix}}.
#' @param x -
#' @param ... -
#' @return A conceptmaps object.
#'@export
conceptmaps <- function(x, ...) UseMethod("conceptmaps")


#' Constructing a conceptmap object
#' 
#' \code{conceptmap} creates an object that encompasses a concept maps. For actual implementations, see
#' \code{\link{conceptmap.default}}, \code{\link{conceptmap.igraph}}, or \code{\link{conceptmap.matrix}}.
#' @param x -
#' @param ... -
#' @return A conceptmap object.
#'@export
conceptmap <- function(x, ...) UseMethod("conceptmap")


#' Modify the concepts of concept maps
#' 
#' \code{modify.concepts} modifies the list of concept of a concept map or of all maps of a set.
#' For actual implementations see \code{\link{modify.concepts.conceptmap}}, or \code{\link{modify.concepts.conceptmaps}}.
#' @param x A conceptmap or conceptmaps object.
#' @param concept.list A list of concepts.
#' @param ... -
#' @return A conceptmaps or conceptmap object.
#'@export
modify.concepts <- function(x, concept.list, ...) UseMethod("modify.concepts")


#' Construct a Pathfinder network from a conceptmap or a concept landscape
#' 
#' \code{pathfinder} creates Pathfinder network. For more information and actual implementations
#'  see \code{\link{pathfinder.matrix}}, \code{\link{pathfinder.conceptmaps}}, or \code{\link{pathfinder.igraph}}.
#' @param data The input data.
#' @param q The q parameter of the Pathfinder algorithm.
#' @param r The r parameter of the Pathfinder algorithm.
#' @param ... -
#' @return The Pathfinder network of the input data.
#'@export
pathfinder <- function(data, q, r, ...) UseMethod("pathfinder")


#' Analyze graph measures of a concept map
#' 
#' \code{analyze.graph.measures} analyzes several graph measures. For actual implementations
#'  see \code{\link{analyze.graph.measures.conceptmap}}, or \code{\link{analyze.graph.measures.igraph}}.
#' @param x A conceptmap.
#' @return A list of several graph measures.
#'@export
analyze.graph.measures <- function(x) UseMethod("analyze.graph.measures")