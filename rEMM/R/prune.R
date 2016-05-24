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


setMethod("prune", signature(x = "EMM"),
	function(x, count_threshold, clusters = TRUE, 
		transitions = FALSE, copy = TRUE, compact = TRUE){

	    if(copy) x <- copy(x)

	    if(clusters && nclusters(x)>0) 
		x <- remove_clusters(x, 
		    rare_clusters(x, count_threshold=count_threshold),
		    copy=FALSE)

	    if(transitions && ntransitions(x)>0) 
		x <- remove_transitions(x, 
		    rare_transitions(x, count_threshold=count_threshold),
		    copy=FALSE)

	    if(compact) x <- compact(x)

	    if(copy) x else invisible(x)
	})
