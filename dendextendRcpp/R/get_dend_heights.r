# Copyright (C) Tal Galili
#
# This file is part of dendextendRcpp.
#
# dendextend is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# dendextend is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#



#' @title Get branches height attr from a dendrogram
#' @export
#' @aliases
#' old_get_branches_heights
#' @description 
#' Get branches height attributes from a dendrogram.
#' 
#' This function is intended to override \code{\link[dendextend]{get_branches_heights}}.
#' Since it is 40-100 times faster.
#' 
#' @param tree a dendrogram object.
#' @param sort logical (TRUE). Should the heights be sorted?
#' @param decreasing logical (FALSE). Should the sort be increasing or 
#' decreasing? Not available for partial sorting
#' @param ... not used.
#' @return A numeric vector of the dendrogram's nodes heights (excluding leaves).
#' @author Tal Galili
#' @seealso \code{\link{labels}}, \code{\link{dendrogram}}, \code{\link{attr}}
#' @examples
#' \dontrun{
#' dend = as.dendrogram(hclust(dist(iris[1:150,-5])))
#' plot(dend)
#' get_height(dend)
#' 
#' 
#' attributes(dend) <- NULL
#' Rcpp_get_dend_heights(dend) # knows to through a warning :)
#' 
#' dend = as.dendrogram(hclust(dist(iris[1:150,-5])))
#' Rcpp_get_dend_heights(dend)
#' get_branches_heights(dend)
#' Rcpp_get_dend_heights(dend,T,F)
#' Rcpp_get_dend_heights(dend,T,F)
#' 
#' # require(dendextend)
#' Rcpp_get_dend_heights(dend)
#' dendextendRcpp_get_branches_heights(dend,sort=F)
#' 
#' \dontrun{
#' require(microbenchmark)
#' microbenchmark(
#'    dendextendRcpp::dendextendRcpp_get_branches_heights(dend),
#'    old_get_branches_heights(dend,sort=F)
#' )
#' # Rcpp is about 40-107 times faster!
#' }
#' }
dendextendRcpp_get_branches_heights <- function (tree, sort = TRUE, decreasing = FALSE, ...) 
{
   height <- Rcpp_get_dend_heights(tree, branches_heights=TRUE, labels_heights= FALSE)
   if (sort) 
      height <- sort(height, decreasing = decreasing)
   return(height)
}




# detach( 'package:dendextendRcpp', unload=TRUE )
# require( 'dendextendRcpp' )
# labels(dend)
