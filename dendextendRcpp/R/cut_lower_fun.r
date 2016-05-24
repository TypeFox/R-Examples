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



#' @title Cut a dendrogram using Rcpp - and run a function on the output
#' @export
#' @aliases
#' old_cut_lower_fun
#' Rcpp_cut_lower
#' @description 
#' Cuts the a tree at height h and returns a list with the FUN function
#' implemented on all the sub trees created by cut at height h.
#' This is used for creating a \link[dendextend]{cutree.dendrogram} function,
#' by using the \code{labels} function as FUN.
#' 
#' This is the Rcpp version of the function, offering a 10-60 times improvement
#' in speed (depending on the tree size it is used on).
#' 
#' @param tree a dendrogram object.
#' @param h a scalar of height to cut the tree by.
#' @param warn logical (FALSE) - should the user be warned if reverting to
#' default? (I set it to FALSE since it can be very noisy sometimes...)
#' @param FUN the function to run (default is "labels")
#' @param ... passed to FUN.
#' @return A list with the output of running FUN on each of the 
#' sub trees derived from cutting "tree"
#' @author Tal Galili
#' @seealso \code{\link{labels}}, \code{\link{dendrogram}},
#' \link[dendextend]{cutree} (in dendextend), \link[stats]{cutree} (in stats)
#' @examples
#' 
#' dend = as.dendrogram(hclust(dist(iris[1:4,-5])))
#' dendextendRcpp::Rcpp_cut_lower(dend, .4)
#' dendextendRcpp::Rcpp_cut_lower(dend, .4, FALSE)
#' # this is really cool!
#' dendextendRcpp_cut_lower_fun(dend, .4, labels)
#' lapply(cut(dend, h = .4)$lower, labels)   
#' dendextendRcpp_cut_lower_fun(dend, .4, order.dendrogram)
#' 
#' \dontrun{
#'    # require(dendextend)
#'    require(dendextendRcpp)
#'    dend_big = as.dendrogram(hclust(dist(iris[1:150,-5])))
#'    require(microbenchmark)
#'    microbenchmark(old_cut_lower_fun(dend_big,.1),
#'                   dendextendRcpp::dendextendRcpp_cut_lower_fun(dend_big,.1),
#'                      times = 100)
#'    # about 7-15 times faster. It is faster the larger the tree is, and the lower h is.
#' }
#' 
dendextendRcpp_cut_lower_fun <- function(tree, h, FUN = labels, warn = FALSE, ...) {
   # is.dendrogram is from dendextend
   if(!dendextend::is.dendrogram(tree)) stop("'tree' needs to be a dendrogram. Aborting the function 'cut_lower_fun'.")
   
   if(is.leaf(tree)) return(list(FUN(tree)))
   # else:
#    dend_and_FUN <- function(x) {
#       class(x) = "dendrogram"
#       FUN(x,...)
#    }
#    return(lapply(Rcpp_cut_lower(tree, h), dend_and_FUN))
   return(lapply(Rcpp_cut_lower(tree, h), FUN))
}




# detach( 'package:dendextendRcpp', unload=TRUE )
# require( 'dendextendRcpp' )
# labels(dend)
