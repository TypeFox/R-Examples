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



#' @title Which height will result in which k for a dendrogram
#' @description
#' This is the Rcpp version of the heights_per_k.dendrogram from the 
#' dendextend R package.
#' It is intended to override that function in order to get a speed gain of
#' ranging from 10 times faster (for a tree with 3 leaves), to 130 times (!) 
#' faster (for a tree with 150 leaves).
#' @export
#' @aliases
#' old_heights_per_k.dendrogram
#' Rcpp_get_dend_heights
#' Rcpp_k_per_height
#' Rcpp_k_per_heights
#' @param tree a dendrogram.
#' @param ... not used.
#' @return a vector of heights, with its names being the k clusters that will
#' result for cutting the dendrogram at each height.
#' 
#' @examples
#' \dontrun{
#' 
#' dend = as.dendrogram(hclust(dist(iris[1:150,-5])))
#' 
#' # plutral heightS
#' dendextendRcpp_heights_per_k.dendrogram(dend,1)
#' dendextendRcpp_heights_per_k.dendrogram(dend,.5)
#' dendextendRcpp_heights_per_k.dendrogram(dend,0.05)
#' dendextendRcpp_heights_per_k.dendrogram(dend,0)
#' dendextendRcpp_heights_per_k.dendrogram(dend,-1)
#' 
#' dendextendRcpp_heights_per_k.dendrogram(dend,c(1,.5,.2,0,-1))
#' 
#' \dontrun{
#' require(microbenchmark)
#' dend = as.dendrogram(hclust(dist(iris[1:150,-5])))
#' dend = as.dendrogram(hclust(dist(iris[1:30,-5])))
#' dend = as.dendrogram(hclust(dist(iris[1:3,-5])))
#' microbenchmark(
#'    #    dendextendRcpp::heights_per_k.dendrogram(dend),
#'    dendextendRcpp::dendextendRcpp_heights_per_k.dendrogram(dend),
#'    dendextend::old_heights_per_k.dendrogram(dend)
#' )
#' # improvment is 10 times faster (in Rcpp) for a tree of size 3
#' # 76 times faster for a tree of size 30
#' # And:
#' # 134 times faster for a tree of size 150!!
#' }
#' }
dendextendRcpp_heights_per_k.dendrogram <- function(tree,...)
{
   # gets a dendro tree
   # returns a vector of heights, and the k clusters we'll get for each of them.
   
   our_dend_heights <- sort(unique(dendextend::get_branches_heights(tree, sort = FALSE)), TRUE)
   # notice that get_branches_heights is a function on dendextendRcpp
   
   heights_to_remove_for_A_cut <- min(-diff(our_dend_heights))/2 # the height to add so to be sure we get a "clear" cut
   heights_to_cut_by <- c((max(our_dend_heights) + heights_to_remove_for_A_cut),   # adding the height for 1 clusters only (this is not mandetory and could be different or removed)
                          (our_dend_heights - heights_to_remove_for_A_cut))
   
   our_dend_ks <- Rcpp_k_per_heights(tree, heights_to_cut_by) # the HUGE improvement over the R code...
   names(heights_to_cut_by) <- our_dend_ks
   return(heights_to_cut_by)
   # notice we might have certion k's that won't exist in this list!
}



# detach( 'package:dendextendRcpp', unload=TRUE )
# require( 'dendextendRcpp' )
# labels(dend)
