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


## build for EMM


## make  newdata a matrix (with a single row)
setMethod("build", signature(x = "EMM", newdata = "numeric"),
	function(x, newdata, verbose = FALSE) build(x, 
		as.matrix(rbind(newdata), verbose))
	)

setMethod("build", signature(x = "EMM", newdata = "data.frame"),
	function(x, newdata, verbose = FALSE) build(x, as.matrix(newdata), 
		verbose)
	)

setMethod("build", signature(x = "EMM", newdata = "matrix"),
	function(x, newdata, verbose = FALSE) {

	    if(verbose) cat("Adding", nrow(newdata) , "observations.","\n")
	    
	    ## cluster all the data (the variable data is in an 
	    ## environment, so there is no need for x <- cluster(x, newdata))
	    cluster(x, newdata, verbose=verbose)

	    ## now update TRACDS (iterate over cluster assignments in last)
	    update(x, last_clustering(x), verbose=verbose)
	    
	    invisible(x)
	}
	)
