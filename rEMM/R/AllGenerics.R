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


## These generics already exist:
## setGeneric("predict", function(object, ...) standardGeneric("predict"))
## setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

## TRACDS
## size is also used in package arules
setGeneric("update", function(object, ...) standardGeneric("update"))
setGeneric("copy", function(x) standardGeneric("copy"))
setGeneric("compact", function(x) standardGeneric("compact"))
setGeneric("size", function(x, ...) standardGeneric("size"))
setGeneric("nclusters", function(x, ...) standardGeneric("nclusters"))
setGeneric("nstates", function(x, ...) standardGeneric("nstates"))
setGeneric("ntransitions", function(x, ...) standardGeneric("ntransitions"))
setGeneric("current_state", function(x) standardGeneric("current_state"))
setGeneric("states", function(x) standardGeneric("states"))
setGeneric("transitions", function(x) standardGeneric("transitions"))
setGeneric("transition", function(x, from, to, ...) 
	standardGeneric("transition"))
setGeneric("transition_matrix", function(x, ...) 
	standardGeneric("transition_matrix"))
setGeneric("transition_table", function(x, newdata, ...) 
	standardGeneric("transition_table"))
setGeneric("initial_transition", function(x, ...) 
	standardGeneric("initial_transition"))
setGeneric("rare_transitions", function(x, count_threshold, ...) 
	standardGeneric("rare_transitions"))
setGeneric("remove_transitions", function(x, from, to, copy=TRUE) 
	standardGeneric("remove_transitions"))
setGeneric("remove_selftransitions", function(x, copy=TRUE) 
	standardGeneric("remove_selftransitions"))
setGeneric("smooth_transitions", function(x, ...) standardGeneric("smooth_transitions"))


## tNN
setGeneric("cluster", function(x, newdata, ...) standardGeneric("cluster"))
setGeneric("clusters", function(x) standardGeneric("clusters"))
setGeneric("cluster_counts", function(x) standardGeneric("cluster_counts"))
setGeneric("last_clustering", function(x,...) standardGeneric("last_clustering"))
setGeneric("cluster_centers", function(x) standardGeneric("cluster_centers"))
setGeneric("find_clusters", function(x, newdata, ...) 
	standardGeneric("find_clusters"))
setGeneric("rare_clusters", function(x, count_threshold, ...) 
	standardGeneric("rare_clusters"))

## EMM
setGeneric("build", function(x, newdata, ...) standardGeneric("build"))
setGeneric("reset", function(x) standardGeneric("reset"))
setGeneric("score", function(x, newdata, ...) standardGeneric("score"))
setGeneric("fade", function(x, t, lambda) standardGeneric("fade"))
setGeneric("prune", function(x, ...) 
	standardGeneric("prune"))
setGeneric("merge_clusters", function(x, to_merge, ...) 
	standardGeneric("merge_clusters"))
setGeneric("remove_clusters", function(x, to_remove, copy=TRUE) 
	standardGeneric("remove_clusters"))

setGeneric("object.size", function(x) standardGeneric("object.size"))

## FIXME: make it one recluster method
setGeneric("recluster_hclust", function(x, ...) 
	standardGeneric("recluster_hclust"))
setGeneric("recluster_kmeans", function(x, ...) 
	standardGeneric("recluster_kmeans"))
setGeneric("recluster_pam", function(x, ...) 
	standardGeneric("recluster_pam"))
setGeneric("recluster_reachability", function(x, ...) 
	standardGeneric("recluster_reachability"))
setGeneric("recluster_tNN", function(x, ...) 
	standardGeneric("recluster_tNN"))
setGeneric("recluster_transitions", function(x, ...) 
	standardGeneric("recluster_transitions"))

