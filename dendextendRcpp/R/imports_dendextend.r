# Copyright (C) Tal Galili
#
# This file is part of dendextend.
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

 
# This file includes functions from "dendextend" that have been modified by 
# this package. But here I keep an "old" copy of them so that they may
# be inspected...


#       assign("old_get_branches_heights", get_branches_heights, envir=as.environment("package:dendextendRcpp"))
#       assign("old_heights_per_k.dendrogram", heights_per_k.dendrogram, envir=as.environment("package:dendextendRcpp"))
#       assign("old_cut_lower_fun", cut_lower_fun, envir=as.environment("package:dendextendRcpp"))


old_cut_lower_fun <- function(tree, h, FUN = labels, warn = FALSE, ...) {
   
   if(!dendextend::is.dendrogram(tree)) stop("'tree' needs to be a dendrogram. Aborting the function 'cut_lower_labels'.")
   
   if(is.leaf(tree)) return(list(FUN(tree)))
   # else:
   #    dend_and_FUN <- function(x) {
   #       class(x) = "dendrogram"
   #       FUN(x,...)
   #    }
   #    return(lapply(Rcpp_cut_lower(tree, h), dend_and_FUN))
   
   return(lapply(cut(tree, h = h)$lower, FUN))   # If the proper labels are not important, this function is around 10 times faster than using labels (so it is much better for some other algorithms)
   
}


old_heights_per_k.dendrogram <- function (tree, ...) 
{
   our_dend_heights <- sort(unique(dendextend::get_branches_heights(tree, 
                                                        sort = FALSE)), TRUE)
   heights_to_remove_for_A_cut <- min(-diff(our_dend_heights))/2
   heights_to_cut_by <- c((max(our_dend_heights) + heights_to_remove_for_A_cut), 
                          (our_dend_heights - heights_to_remove_for_A_cut))
   names(heights_to_cut_by) <- sapply(heights_to_cut_by, function(h) {
      length(cut(tree, h = h)$lower)
   })
   names(heights_to_cut_by)[1] <- "1"
   return(heights_to_cut_by)
}



old_get_branches_heights <- function (tree, sort = TRUE, decreasing = FALSE, ...) 
{
   height <- dendextend::get_nodes_attr(tree, "height", include_leaves = FALSE, 
                            na.rm = TRUE)
   if (sort) 
      height <- sort(height, decreasing = decreasing)
   return(height)
}


