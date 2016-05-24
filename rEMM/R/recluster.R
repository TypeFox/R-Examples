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


## FIXME: find medoids for hclust

## medoid is defined as the object of a cluster, whose average 
## dissimilarity to all the objects in the cluster is minimal
## min_{m\inC}(1/n_C sum_{i\inC\m}(d(i,m))
.find_medoids <- function(d, k, cl) {
    dm <- as.matrix(d)
    sapply(1:k, FUN =function(i){
		take <- cl==i    
		names(which.min(colSums(dm[take,take, drop=FALSE])))
	    })
}

## hierarchical clustering
setMethod("recluster_hclust", signature(x = "EMM"),
	function(x, k=NULL, h=NULL,  method="average", 
		..., prune=NULL, copy=TRUE) {

	    if(copy) x <- copy(x)

	    if(!is.null(prune)) x <- prune(x, count_threshold = prune, 
		    transitions = FALSE, copy = FALSE)

	    d <- dist(cluster_centers(x), method = x@distFun)
	    hc <- hclust(d, method=method, ...)
	    cl <- cutree(hc, k=k, h=h)

	    ## if only h was given
	    if(is.null(k)) k <- max(cl)

	    if(is(cl, "matrix")){ 
		x <- lapply(1:ncol(cl), FUN=function(i){
			    if(!x@centroids) 
				new_center <- cluster_centers(x)[.find_medoids(d, k, cl[,i]),]
			    ## centroids are handled by merge_clusters!
			    else 
				new_center <- NULL
				
			    merge_clusters(x, cl[,i], 
				    clustering=TRUE, 
				    new_center = new_center, 
				    copy=TRUE)
			})
	    }else{ 
		if(!x@centroids) 
		    new_center <- cluster_centers(x)[.find_medoids(d, k, cl),]
		else 
		    new_center <- NULL
		    
		merge_clusters(x, cl, 
			clustering=TRUE,  
			new_center = new_center,
			copy=FALSE)
	    }

	    attr(x, "cluster_info") <- list(clustering=cl, dendrogram=hc)
	    
	    if(copy) x
	    else invisible(x)
	}
	)

## k-means (euclidean)
setMethod("recluster_kmeans", signature(x = "EMM"),
	function(x, k, ..., prune=NULL, copy=TRUE) {
	    
	    if(copy) x <- copy(x)

	    if(!is.null(prune)) x <- prune(x, count_threshold = prune, 
		    transitions = FALSE, copy = FALSE)

	    if(!identical(tolower(x@measure), "euclidean")) warning(
		    paste("Using k-means implies Euclidean distances but the EMM uses:", 
			    x@measure))

	    cl <- kmeans(cluster_centers(x), centers = k, ...)

	    merge_clusters(x, cl$cluster, 
		    clustering=TRUE, 
		    new_center=cl$centers,
		    copy=FALSE)

	    attr(x, "cluster_info") <- cl
	    
	    if(copy) x
	    else invisible(x)
	}
	)

## Partitioning around medoids (k-medians)
setMethod("recluster_pam", signature(x = "EMM"),
	function(x, k, ..., prune=NULL, copy=TRUE) {

	    if(copy) x <- copy(x)
	    if(!is.null(prune)) x <- prune(x, count_threshold = prune, 
		    transitions = FALSE, copy = FALSE)

	    d <- dist(cluster_centers(x), method = x@distFun)
	    cl <- pam(d, k=k, ...)

	    medoids <- cluster_centers(x)[cl$medoids,]
	    
	    merge_clusters(x, cl$clustering, 
		    clustering=TRUE, 
		    new_center=medoids,
		    copy=FALSE)

	    attr(x, "cluster_info") <- cl
	    
	    if(copy) x
	    else invisible(x)
	}
	)

## reachability
setMethod("recluster_reachability", signature(x = "EMM"),
	function(x, h, ..., prune=NULL, copy=TRUE) {


	    if(copy) x <- copy(x)
	    
	    if(!is.null(prune)) x <- prune(x, count_threshold = prune, 
		    transitions = FALSE, copy = FALSE)

	    d <- as.matrix(dist(cluster_centers(x), method = x@distFun))

	    # get adjecency matrix and find all paths 
	    a_mat <- d < h
	    r_mat <- a_mat
	    for(i in 1:size(x)) {
		r_mat <- r_mat%*%a_mat
		storage.mode(r_mat) <- "logical"
	    }

	    to_merge <- unique(apply(r_mat, MARGIN=1, FUN = which))
	    to_merge <- lapply(to_merge, as.character)

	    for(i in 1:length(to_merge)) {
		m <- to_merge[[i]]
		if(length(m)>1) {
		    merge_clusters(x, to_merge = m, copy=FALSE)
		}
	    }

	    if(copy) x
	    else invisible(x)
	})


## tNN
setMethod("recluster_tNN", signature(x = "EMM"),
	function(x, threshold=NULL, ..., prune=NULL, copy=TRUE) {
	    
	    if(copy) x <- copy(x)

	    if(!is.null(prune)) x <- prune(x, count_threshold = prune, 
		    transitions = FALSE, copy = FALSE)
	    
	    if(is.null(threshold)) threshold <- x@threshold

	    cl <- tNN(threshold=threshold, measure=x@measure, 
		    centroids=x@centroids, lambda=0) 

	    cluster(cl, cluster_centers(x))
	    assignments <- last_clustering(cl)
	    
	    merge_clusters(x, as.integer(assignments), clustering=TRUE, 
		    copy=FALSE)

	    if(copy) x
	    else invisible(x)
	})

## transitions: group all states which intersecting radius
## and then homogenizes transitions between groups. 
## Note: does not cluster states!
setMethod("recluster_transitions", signature(x = "EMM"),
	function(x, threshold=NULL, ..., prune=NULL, copy=TRUE) {
	   
	    if(is.null(threshold)) threshold <- 2*x@threshold

	    if(copy) x <- copy(x)

	    if(!is.null(prune)) x <- prune(x, count_threshold = prune, 
		    transitions = FALSE, copy = FALSE)

	    d <- dist(cluster_centers(x), method = x@distFun)
	    hc <- hclust(d, method="single")
	    
	    ## FIXME: this ignores vat_thresholds for now!
	    clusters <- cutree(hc, h=threshold)

	    for(i in unique(clusters)) {
		cl <- which(clusters==i)
		ncl <- which(clusters!=i)

		if(length(cl)>1) {
		    ### fix transitions between elements in group
		    x@tracds_d$mm@counts[cl,cl] <- sum(
			    x@tracds_d$mm@counts[cl,cl])/length(cl)	

		    ### fix outgoing transitions
		    x@tracds_d$mm@counts[cl,ncl] <- matrix(
			    colSums(x@tracds_d$mm@counts[cl,ncl])/length(cl), 
			    nrow=length(cl), ncol=length(ncl), byrow=TRUE)
		    
		    ### fix incoming transitions
		    x@tracds_d$mm@counts[ncl,cl] <- matrix(
			    rowSums(x@tracds_d$mm@counts[ncl,cl])/length(cl), 
			    nrow=length(ncl), ncol=length(cl), byrow=FALSE)

		}
	    }


	    if(copy) x
	    else invisible(x)
	})

