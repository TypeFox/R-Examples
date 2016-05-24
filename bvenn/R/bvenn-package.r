#' A Simple alternative to proportional Venn diagrams
#'
#' A Venn diagram is a useful tool to visualise the overlap between sets of elements. 
#' However, often the proportions of the overlaps are not accurately depicted by the 
#' areas on the figure. In fact, for some configurations of the 3 sets it is even 
#' impossible. This means that the proportions of the sets and overlaps often still have 
#' to be deduced from the numbers associated with the figure. 
#' 
#' This package implements a simple alternative to the traditional Venn diagram, where we 
#' depict each overlap as a separate bubble with area proportional to the overlap size. 
#' Relation of the bubbles to input sets becomes clear from their arrangement.  
#' 
#' These figures are much easier to read than normal approximately proportional Venn 
#' diagrams. First, the sets are always in the same positions on the plot and not moving 
#' around like on normal Venn diagrams. This feature makes it easy to compare the figures 
#' if several diagrams are drawn side-by-side. Second, the proportions are easier to 
#' read, because they all are depicted as similar shapes not as circles and cuts of 
#' circles as on ordinary Venn diagrams.
#' 
#' The function for drawing the figures is \code{\link{bvenn}}.
#' 
#' @name bvenn-package
#' @docType package
#' 
#' @import grid
NA