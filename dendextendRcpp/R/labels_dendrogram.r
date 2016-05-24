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



#' @title Find Labels from a dendrogram Object using Rcpp
#' @export
#' @aliases 
#' Rcpp_labels_dendrogram
#' stats_labels.dendrogram
#' @description 
#' Extract the leaves labels from a dendrogram object.
#' @param object a dendrogram object.
#' @param warn logical (FALSE) - should the user be warned if reverting to
#' default? (I set it to FALSE since it can be very noisy sometimes...)
#' @param ... not used.
#' @return A vector of labels from the dendrogram leaves.
#' This is often a character vector, but there are cases it might be integer.
#' @author Romain Francois, Dirk Eddelbuettel, Tal Galili
#' @source 
#' R-devel-mailing list.
#' @seealso \code{\link{labels}}, \code{\link{dendrogram}}
#' @examples
#' dend <- as.dendrogram(hclust(dist(USArrests)))
#' 
#' labels(dend)
#' 
#' \dontrun{
#' # require(microbenchmark)
#' microbenchmark::microbenchmark(dendextendRcpp::stats_labels.dendrogram(dend),
#'                dendextendRcpp::dendextendRcpp_labels.dendrogram(dend),
#'                times = 100)
#' # about 30 times faster. It is faster the larger the tree is.
#' }
dendextendRcpp_labels.dendrogram <- function(object, warn = FALSE, ...) {
   
   if(is.leaf(object)) return(attr(object, "label"))   
   
   # we would get errors if, for example, labels are not characters
   tryCatch(return(Rcpp_labels_dendrogram(object)) , error = function(e) {
      if(warn) warning("Your tree's labels are not 'character'. Hence the 'labels' function can not use the Rcpp function \n and is expected to be 20 to 40 times SLOWER! \n In order to fix this, simply run on your tree:\n labels(tree)<-as.character(labels(tree)) \n This function is in the {dendextend} R package \n Do it once - and many functions which rely on the tree 'labels' will run faster. :) ")
   })
   # else: (in case of integer labels)
   return(stats_labels.dendrogram(object))
}


# labels.dendrogram <- dendextendRcpp_labels.dendrogram


# detach( 'package:dendextendRcpp', unload=TRUE )
# require( 'dendextendRcpp' )
# labels(dend)
