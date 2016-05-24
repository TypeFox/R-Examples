#######################################################################
# rEMM - Extensible Markov Model (EMM) for Data Stream Clustering in R
# Copyrigth (C) 2011 Michael Hahsler
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


## clustering = TRUE gets integer in to_merge
setMethod("merge_clusters", signature(x = "EMM", to_merge = "integer"),
	function(x, to_merge, clustering = FALSE, new_center = NULL, 
		copy=TRUE) {

	    if(copy) x <- copy(x)

	    ## handle a clustering
	    if(!clustering) stop("to_merge needs to be all character!")
	    k <- max(to_merge)

	    if(!is.null(new_center) && nrow(new_center) != k) 
		stop("new_center has not the right number of rows.")

	    orig_states <- clusters(x)
	    
	    for(i in 1:k) {
		m <- orig_states[to_merge==i]
		if(length(m)>1)
		    merge_clusters(x, m, clustering = FALSE, 
			new_center[i,], copy = FALSE)
	    }

	    invisible(x)
	})

## clustering = FALSE gets character
setMethod("merge_clusters", signature(x = "EMM", to_merge = "character"),
	function(x, to_merge, clustering = FALSE, new_center = NULL, 
		copy = TRUE) {

	    if(copy) x <- copy(x)
	    
	    if(clustering) 
		stop("to_merge has the wrong format for clustering!")

	    if(!all(is.element(to_merge, clusters(x)))) 
		stop("not all clusters in to_merge exist in x!")
	    
	    ## nothing to do
	    if(length(to_merge) < 2) return(x)

	    new_state <- to_merge[1]
	    to_delete <- states(x) %in% to_merge[-1]


	    ## TRACDS
	    x@tracds_d$mm <- smc_mergeStates(x@tracds_d$mm, to_merge)

	    if(x@tracds_d$current_state %in% to_delete) 
		x@tracds_d$current_state <- new_state

	    ## tNN
	    ## save old state centers
	    old_centers <- cluster_centers(x)[to_merge,]

	    ## create new state
	    if(is.null(new_center)) {
		if(x@centroids) {
		    x@tnn_d$centers[new_state,] <- 
		    colSums(old_centers*x@tnn_d$counts[to_merge])/
		    sum(x@tnn_d$counts[to_merge])
		}else {
		    ## we take the medoid of the larger cluster
		    x@tnn_d$centers[new_state,] <-
		    old_centers[which.max(cluster_counts(x)[to_merge]),]
		}
	    }else{ 
		## user supplied new center
		if(identical(length(new_center), ncol(x@tnn_d$centers))) 
		    x@tnn_d$centers[new_state,] <- new_center
		else stop("new_center does not have the correct length/ncol!")
		}


	    x@tnn_d$counts[new_state] <- sum(x@tnn_d$counts[to_merge])

	    x@tnn_d$centers <- x@tnn_d$centers[!to_delete, , drop=FALSE]
	    x@tnn_d$counts <- x@tnn_d$counts[!to_delete]



	    ## FIXME: this only works for metric dissimilarities (distances)
	    ## new threshold is max. dissimilarity vom new centroid to any old
	    ## centroid + its threshold

	    d <- dist(cluster_centers(x)[new_state,,drop=FALSE], 
		    old_centers, method=x@measure)[1,]

	    new_threshold <- max(d + x@tnn_d$var_thresholds[names(d)])
	    names(new_threshold) <- new_state

	    x@tnn_d$var_thresholds[new_state] <- new_threshold

	    ## remove var. thresholds
	    x@tnn_d$var_thresholds <- x@tnn_d$var_thresholds[!to_delete]

	    invisible(x)
	})

