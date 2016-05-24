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


## create a tNN clustering from k-means, etc.
## This is TRAC (without DS)

TRAC <- function(x, data = NULL, measure = "euclidean") {
    if(is(x, "kmeans")) {
	counts <- x$size
	k <- length(counts)
	centers <- x$centers
	thresholds <- sqrt(x$withinss/x$size)*3
	order <- x$cluster
    
    ### PAM
    }else if(is(x, "partition")) {
	counts <- x$clusinfo[,"size"]
	k <- length(counts)
	centers <- x$medoids
	thresholds <- x$clusinfo[,"max_diss"]
	order <- x$cluster
    
    ### x must be a cluster index vector
    }else if(is.numeric(x) && !is.null(data)) {
	order <- as.integer(x)
	k <- max(order)
	counts <- tabulate(order, k)
	centers <- t(sapply(1:k, 
			FUN=function(i) colMeans(data[order==i,,drop=FALSE])))
	thresholds <- sapply(1:k,
		        FUN=function(i) max(dist(centers[i,,drop=FALSE], 
				data[order==i,,drop=FALSE], method = measure)))

    }else
	stop("Needs to be a kmeans or partition (PAM) object or a cluster vector and the used data.")

    emm <- new("EMM", measure=measure, threshold = max(thresholds))

    ## create tNN
    states <- as.character(1:k)
    names(counts) <- states

    emm@tnn_d$centers <- centers
    emm@tnn_d$counts <- counts
    emm@tnn_d$var_thresholds <- thresholds

    ## update TRACDS
    ## make sure the order of states corresponds to tNN
    emm@tracds_d$mm <- smc_addState(emm@tracds_d$mm, states)
    update(emm, order)

    emm
}

