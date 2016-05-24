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


setMethod("remove_clusters", signature(x = "EMM", to_remove = "character"),
	function(x, to_remove, copy=TRUE) {

		if(copy) x <- copy(x)

		if(length(to_remove)==0) return(x)
		
		to_remove_pos <- states(x) %in% to_remove

		## TRACDS 
		x@tracds_d$mm <- smc_removeState(x@tracds_d$mm, to_remove)
		if(is.element(x@tracds_d$current_state, to_remove)) 
		    x@tracds_d$current_state <- as.character(NA)

		## tNN
		x@tnn_d$centers <- x@tnn_d$centers[!to_remove_pos,]
		x@tnn_d$counts <- x@tnn_d$counts[!to_remove_pos]
		x@tnn_d$var_thresholds <- x@tnn_d$var_thresholds[!to_remove_pos]

		invisible(x)
	}
)

setMethod("remove_transitions", signature(x = "EMM", 
		from ="matrix", to="missing"),
	function(x, from, to, copy=TRUE) 
	remove_transitions(x, from[,1], from[,2], copy)
)

setMethod("remove_transitions", signature(x = "EMM", 
		from ="character", to="character"),
	function(x, from, to, copy=TRUE) {

	    if(copy) x <- copy(x)

	    if(length(from) != length(to)) 
		stop("length of from and to do not match!")
	    if(length(from)==0) return(x)

	    x@tracds_d$mm <- smc_removeTransition(x@tracds_d$mm,from, to)
	    
	    if(copy) x else invisible(x)
	}
	)

setMethod("remove_selftransitions", signature(x = "EMM"),
	function(x, copy=TRUE) {

	    if(copy) x <- copy(x)

	    x@tracds_d$mm <- smc_removeSelfTransition(x@tracds_d$mm)	
	    
	    if(copy) x else invisible(x)
	}
	)
